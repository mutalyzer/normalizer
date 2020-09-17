import json
import time
from functools import lru_cache

import extractor
from mutalyzer_hgvs_parser import parse_description_to_model
from mutalyzer_mutator.mutator import mutate
from mutalyzer_retriever import retriever

from .converter import (
    de_to_hgvs,
    to_cds_coordinate,
    to_delins,
    to_hgvs_locations,
    to_internal_locations,
)
from .converter.to_internal_coordinates import location_to_internal_coordinate
from .converter.to_internal_indexing import to_internal_indexing
from .description import get_coordinate_system, get_selector_id, model_to_string, variant_to_description
from .protein import get_protein_description, get_protein_descriptions
from .reference import (
    extract_sequences,
    get_available_selectors,
    get_mol_type,
    get_selector_model,
    get_selectors_ids,
    get_reference_model
)
from .util import get_time_information, string_k_v
from .visualization import to_be_visualized
from .checker import run_checks


def extract_sequences(references):
    """
    Return a dictionary with reference ids as keys and their corresponding
    sequences as values.

    :param references: Dictionary with reference models.
    :rtype: dict
    :return: Reference ids as keys and their corresponding sequences as values
    """
    sequences = {}
    for reference in references:
        sequences[reference] = references[reference]["sequence"]["seq"]
    return sequences


class Description(object):
    def __init__(self, description):
        self.input_description = description
        self._input_description_model = None
        self._normalized_description = None
        self._equivalent_descriptions = []
        self._reference_id = None
        self._reference_models = {}
        self._selector_id = None
        self._selector_model = None
        self._coordinate_system = None
        self._mol_type = None
        self._de_hgvs_variants = None
        self.status = {"errors": [], "warnings": []}
        self._time_stamps = []

        self._normalize()

    def __repr__(self):
        output = "{}\n".format(self.input_description)
        w = 20
        if self._input_description_model:
            output += string_k_v(w, "Syntax check", "Pass")
        else:
            output += string_k_v(w, "Syntax check", "Failed")
        if self._reference_id:
            output += string_k_v(w, "Reference ID", self._reference_id)
        if self._coordinate_system:
            output += string_k_v(w, "Coordinate system", self._coordinate_system)
        if self._reference_models:
            output += string_k_v(w, "Reference model", "Retrieved")
        else:
            output += string_k_v(w, "Reference model", "Not retrieved")
        if self._mol_type:
            output += string_k_v(w, "Reference mol type", self._mol_type)
        else:
            output += string_k_v(w, "Reference mol type", "Not retrieved")
        if self._selector_id:
            output += string_k_v(w, "Selector ID", self._selector_id)
            if self._selector_model:
                output += string_k_v(w, "Selector model", "Retrieved")
            else:
                output += string_k_v(w, "Selector model", "Not retrieved")
        else:
            output += string_k_v(w, "Selector ID", "-")
        if self.status["errors"]:
            output += " Errors:\n"
            for error in self.status["errors"]:
                output += "  - {}\n".format(error)
        if self.status["warnings"]:
            output += " Warnings:\n"
            for warning in self.status["warnings"]:
                output += "  - {}\n".format(warning)
        return output

    def _add_error(self, error):
        self.status["errors"].append(error)

    def _add_warning(self, warning):
        self.status["warnings"].append(warning)

    def syntax_check(self):
        """
        Calls the HGVS parser and retrieves the description model.
        If successful, the self._input_description_model is populated,
        otherwise the parsing error is added to the errors list.
        """
        try:
            model = parse_description_to_model(self.input_description)
        except Exception as e:
            model = {"errors": "Some parsing error occured."}
        if model.get("errors"):
            self._add_error(model["errors"])
        else:
            self._input_description_model = model

    def construct_reference(self):
        """
        Populates the instance reference attributes.
        The following steps are performed:
        - gets the reference ID from the description model.
        - retrieves the reference models.
        - retrieves the selector ID and its feature mode.
        - identifies the reference molecule type.
        """
        self._reference_id = self._input_description_model["reference"]["id"]
        self._append_reference(self._reference_id)
        self._reference_models["reference"] = self._reference_models[self._reference_id]
        self._selector_id = get_selector_id(self._input_description_model)
        if self._selector_id:
            self._selector_model = get_selector_model(
                self._reference_models[self._reference_id]["model"], self._selector_id
            )
        self._coordinate_system = get_coordinate_system(self._input_description_model)
        self._mol_type = get_mol_type(self._reference_models[self._reference_id])

    def _handle_no_description_selector(self):
        available_selectors = get_available_selectors(
            self._reference_models[self._reference_id]["model"], self._coordinate_system
        )
        if len(available_selectors) == 0:
            self._add_error(
                "ENOSELECTOR: {} coordinate system used but no selector ID "
                "provided in the description. In addition, there is no "
                "selector available in the reference model.".format(
                    self._coordinate_system
                )
            )
        elif len(available_selectors) == 1:
            self._add_warning(
                "WNOSELECTOR: {} coordinate system used but no selector ID "
                "provided in the description. Only {} present in the reference,"
                " which is chosen as default.".format(
                    self._coordinate_system, available_selectors[0]
                )
            )
            self._selector_id = available_selectors[0]
            self._selector_model = get_selector_model(
                self._reference_models[self._reference_id]["model"], self._selector_id
            )
        elif len(available_selectors) > 1:
            self._add_error(
                "ENOSELECTOR: {} coordinate system used but no selector ID "
                "provided in the description. Please choose between the "
                "following selectors available in the reference: {}".format(
                    self._coordinate_system, available_selectors
                )
            )

    def check_description_reference_consistency(self):
        if self._mol_type in ["dna", "genomic DNA"]:
            if self._coordinate_system in ["c", "n"]:
                if self._selector_id is None:
                    self._handle_no_description_selector()

    def _append_reference(self, reference_id):
        """
        Retrieves and appends to the _reference_models the reference model
        corresponding to provided reference_id.
        """
        if reference_id not in self._reference_models.keys():
            reference = get_reference_model(reference_id)
            if reference is None:
                self._add_error(
                    "No reference was retrieved for {}.".format(reference_id)
                )
            else:
                self._reference_models[reference_id] = reference

    def _get_sequence(self, reference_id):
        """
        Retrieves from the _reference_models the sequence that corresponds to
        the provided reference_id.
        :param reference_id:
        :return:
        """
        return self._reference_models[reference_id]["sequence"]["seq"]

    def _append_reference_sequence(self, reference_id, sequence):
        """
        Appends the provided sequence in the _reference_models.
        """
        self._reference_models[reference_id] = {"sequence": {"seq": sequence}}

    def _get_sequences(self):
        """
        Retrieves a dictionary from the _reference_models with reference ids
        as keys and their corresponding sequences as values.
        """
        return extract_sequences(self._reference_models)

    def _mutate(self):

        self._append_reference_sequence(
            "observed", mutate(self._get_sequences(), self._delins_variants)
        )

    def get_equivalent_descriptions(self):
        equivalent_descriptions = []

        transcript_ids = get_selectors_ids(
            self._reference_models[self._reference_id]["model"]
        )

        for transcript_id in transcript_ids:
            equivalent_variant_model = to_hgvs_locations(
                self._de_hgvs_variants,
                self._reference_models[self._reference_id],
                transcript_id,
            )

            equivalent_descriptions.append(model_to_string(equivalent_variant_model))
        return equivalent_descriptions

    def _get_normalized_description(self):
        self._normalized_description_model = to_hgvs_locations(
            self._de_hgvs_variants,
            self._reference_models[self._reference_id],
            self._selector_id,
            degenerate=True
        )

        self.normalized_description = model_to_string(
            self._normalized_description_model
        )

    def _is_to_extract_protein_description(self):
        if self._coordinate_system == "c" and self._selector_model:
            return True

    def _normalize(self):
        def print_vars(vars):
            v_s = []
            for v in vars:
                v_s.append(variant_to_description(v))
            print('; '.join(v_s))

        self._time_stamps.append(("initial", time.perf_counter()))

        self.syntax_check()
        print(json.dumps(self._input_description_model, indent=2))

        if self.status["errors"]:
            return

        self._time_stamps.append(("syntax parser", time.perf_counter()))

        self.construct_reference()

        self.check_description_reference_consistency()

        self._time_stamps.append(("retriever", time.perf_counter()))

        internal_locations_variants = to_internal_locations(
            self._input_description_model, self._reference_models)["variants"]

        stop, messages = run_checks(internal_locations_variants, self._reference_models)
        self.status["messages"] = messages
        if stop:
            return

        if self.status["errors"]:
            return
        self._delins_variants = to_delins(internal_locations_variants)

        internal_coordinate = ''

        print_vars(self._delins_variants)

        self._time_stamps.append(("to delins", time.perf_counter()))

        self._mutate()
        self._time_stamps.append(("mutator", time.perf_counter()))

        de_variants = extractor.describe_dna(
            self._get_sequence("reference"), self._get_sequence("observed")
        )

        self._time_stamps.append(("description extractor", time.perf_counter()))

        self._de_hgvs_variants = de_to_hgvs(de_variants, self._get_sequences())

        self._get_normalized_description()

        # self.status["equivalent descriptions"] = self.get_equivalent_descriptions()

        self.status["normalized description"] = self.normalized_description

        self._time_stamps.append(("to HGVS description", time.perf_counter()))

        self.status["time information (s)"] = get_time_information(self._time_stamps)

        # self.status["visualize"] = to_be_visualized(
        #     de_variants, self._normalized_description_model["variants"]
        # )

        # protein_descriptions = get_protein_descriptions(
        #     de_variants, self._reference_models
        # )
        # if protein_descriptions:
        #     self.status["protein descriptions"] = protein_descriptions


def mutalyzer3(hgvs_description):

    description = Description(hgvs_description)

    return {
        k: description.status[k] for k in description.status if description.status[k]
    }
