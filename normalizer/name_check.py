import copy

from extractor import describe_dna
from mutalyzer_mutator import mutate

from .converter.to_delins import to_delins
from .converter.to_hgvs import to_hgvs_locations
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description import (
    description_to_model,
    get_errors,
    get_references_from_description_model,
    model_to_string,
)


class Description(object):
    def __init__(self, description):
        self.input_description = description
        self.augmented_description = None
        self.normalized_description = None
        self.input_model = description_to_model(description)
        self.augmented_model = {}
        self.internal_coordinates_model = {}
        self.internal_indexing_model = {}
        self.delins_model = {}
        self.de_model = {}
        self.de_hgvs_internal_indexing_model = {}
        self.de_hgvs_coordinate_model = {}
        self.de_hgvs_model = {}
        self.references = {}
        self.observed_sequence = None

    def augment_input_model(self):
        self.augmented_model = copy.deepcopy(self.input_model)
        if get_errors(self.augmented_model):
            return
        get_references_from_description_model(self.augmented_model, self.references)
        if get_errors(self.augmented_model):
            return
        else:
            self.augmented_description = model_to_string(self.augmented_model)

    def get_internal_coordinate_model(self):
        if self.augmented_model and not get_errors(self.augmented_model):
            self.internal_coordinates_model = to_internal_coordinates(
                self.augmented_model
            )

    def get_internal_indexing_model(self):
        if self.internal_coordinates_model and not get_errors(
            self.internal_coordinates_model
        ):
            self.internal_indexing_model = to_internal_indexing(
                self.internal_coordinates_model
            )

    def get_delins_model(self):
        if self.internal_indexing_model and not get_errors(
            self.internal_indexing_model
        ):
            self.delins_model = to_delins(self.internal_indexing_model)

    def _get_sequences(self):
        """
        Retrieves a dictionary from the _reference_models with reference ids
        as keys and their corresponding sequences as values.
        """
        sequences = {k: self.references[k].sequence() for k in self.references}
        sequences["reference"] = self.references[
            self.augmented_model["reference"]["id"]
        ].sequence()
        return sequences

    def _mutate(self):
        if self.delins_model and not get_errors(self.delins_model):
            self.observed_sequence = mutate(
                self._get_sequences(), self.delins_model["variants"]
            )

    def extract(self):
        reference_sequence = self.references[
            self.augmented_model["reference"]["id"]
        ].sequence()
        if self.observed_sequence and reference_sequence:
            de_variants = describe_dna(reference_sequence, self.observed_sequence)
            if de_variants:
                self.de_model = {
                    "reference": copy.deepcopy(
                        self.internal_indexing_model["reference"]
                    ),
                    "coordinate_system": "i",
                    "variants": de_variants,
                }

    def get_de_hgvs_internal_indexing_model(self):
        if self.de_model:
            self.de_hgvs_internal_indexing_model = {
                "reference": copy.deepcopy(self.internal_indexing_model["reference"]),
                "coordinate_system": "i",
                "variants": de_to_hgvs(
                    self.de_model["variants"],
                    {
                        "reference": self.references[
                            self.augmented_model["reference"]["id"]
                        ].sequence(),
                        "observed": self.observed_sequence,
                    },
                ),
            }

    def get_de_hgvs_coordinates_model(self):
        if self.augmented_model["reference"].get("selector"):
            selector_id = self.augmented_model["reference"]["selector"]
        else:
            selector_id = None
        if self.de_hgvs_internal_indexing_model:
            self.de_hgvs_model = to_hgvs_locations(
                self.de_hgvs_internal_indexing_model["variants"],
                self.references[self.augmented_model["reference"]["id"]].model,
                selector_id,
                True,
            )

    def get_normalized_description(self):
        if self.de_hgvs_model:
            self.normalized_description = model_to_string(self.de_hgvs_model)

    def normalize(self):
        self.augment_input_model()
        self.get_internal_coordinate_model()
        self.get_internal_indexing_model()
        self.get_delins_model()
        self._mutate()
        self.extract()
        self.get_de_hgvs_internal_indexing_model()
        self.get_de_hgvs_coordinates_model()
        self.get_normalized_description()

    def output(self):
        output = {
            "input_model": self.input_model,
            "augmented_model": self.augmented_model,
            "internal_coordinates_model": self.internal_coordinates_model,
            "internal_indexing_model": self.internal_indexing_model,
            "reference_ids": list(self.references.keys()),
        }
        if self.augmented_description:
            output["augmented_description"] = self.augmented_description
        if self.normalized_description:
            output["normalized_description"] = self.normalized_description
        return output


def normalize(description_to_normalize):
    description = Description(description_to_normalize)
    description.normalize()
    return description.output()