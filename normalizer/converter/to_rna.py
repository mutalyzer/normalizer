import bisect
from copy import deepcopy

from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_mutator.util import reverse_complement

from normalizer.util import get_end, get_start

from ..description_model import variant_to_description
from ..reference import (
    extract_feature_model,
    get_feature,
    get_selector_model,
    slice_to_selector,
)
from .to_hgvs_coordinates import genomic_to_point, reverse_strand_shift


def _yield_locations(annotations):
    yield annotations["location"], annotations["type"]
    if annotations.get("features"):
        for feature in annotations["features"]:
            yield from _yield_locations(feature)


def to_rna_reference_model(reference_model, selector_id):
    """
    Get the RNA reference model of the provided selector.

    1. Extract the tree corresponding to the selector from the model (including
    the parents).
    2. Slice the sequence.
    3. Update the model features locations using the crossmapper.

    TODO: Convert sequence to RNA?
    TODO: Make sure everything is on the plus strand?

    :param reference_model: Reference model.
    :param selector_id: Selector ID.
    :return: RNA reference model.
    """
    rna_model = {
        "annotations": deepcopy(
            extract_feature_model(reference_model["annotations"], selector_id)[0]
        ),
        "sequence": {"seq": slice_to_selector(reference_model, selector_id)},
    }
    s_m = get_selector_model(rna_model["annotations"], selector_id, True)
    x = NonCoding(s_m["exon"]).coordinate_to_noncoding

    new_start = x(s_m["exon"][0][0])[0] - 1
    new_end = x(s_m["exon"][-1][-1])[0]
    for location, f_type in _yield_locations(rna_model["annotations"]):
        if f_type in ["exon", "CDS"]:
            location["start"]["position"] = x(location["start"]["position"])[0] - 1
            location["end"]["position"] = x(location["end"]["position"])[0]
        else:
            location["start"]["position"] = new_start
            location["end"]["position"] = new_end

    return rna_model


def _get_location_type(location, exons):
    start_i = bisect.bisect_right(exons, get_start(location))
    end_i = bisect.bisect_left(exons, get_end(location))
    print(start_i, end_i)
    if abs(start_i - end_i) == 0:
        if start_i % 2 == 0:
            return "same intron"
        else:
            return "same exon"
    elif abs(start_i - end_i) == 1:
        if start_i % 2 == 0:
            return "adjacent intron exon"
        else:
            return "adjacent exon intron"
    elif abs(start_i - end_i) > 1:
        if start_i % 2 == 0 and end_i % 2 == 0:
            return "intron intron"
        elif start_i % 2 == 0 and end_i % 2 == 1:
            return "intron exon"
        elif start_i % 2 == 1 and end_i % 2 == 0:
            return "exon intron"
        elif start_i % 2 == 1 and end_i % 2 == 1:
            return "exon exon"


def to_rna_variants(variants, sequences, selector_model):
    """
    Convert coordinate delins variants to RNA.

    :param variants: Variants with coordinate locations.
    :param sequences: Sequences dictionary.
    :param selector_model: Selector model.
    """
    exons = [e for l in selector_model["exon"] for e in l]
    print(exons)
    for variant in variants:
        print(f"- {variant_to_description(variant)}")
        print(f"  {_get_location_type(variant['location'], exons)}")


def _point_to_cds_coordinate(point, selector_model, crossmap):
    genomic_to_coordinate = Genomic().genomic_to_coordinate
    if selector_model.get("inverted"):
        if point.get("shift"):
            point["position"] -= point["shift"]
    coding = crossmap.coordinate_to_coding(point["position"], degenerate=True)
    if coding[2] == -1:
        return genomic_to_point(0)
    else:
        return genomic_to_point(genomic_to_coordinate(coding[0]))


def get_inserted_sequence(insertion, sequences):
    if isinstance(insertion["source"], str):
        source = insertion["source"]
    elif isinstance(insertion["source"], dict):
        source = insertion["source"]["id"]
    return sequences[source][
        get_start(insertion["location"]) : get_end(insertion["location"])
    ]


def merge_inserted_to_string(inserted, sequences):
    inserted_value = ""
    for insertion in inserted:
        if insertion.get("sequence"):
            inserted_value += insertion.get("sequence")
        else:
            inserted_value += get_inserted_sequence(insertion, sequences)
        if insertion.get("inverted"):
            inserted_value = reverse_complement(inserted_value)

    return {"source": "description", "sequence": inserted_value}


def variant_to_cds_coordinate(variant, sequences, selector_model, crossmap):
    new_variant = deepcopy(variant)

    location = new_variant["location"]

    if location["type"] == "range":
        location["start"] = _point_to_cds_coordinate(
            location["start"], selector_model, crossmap
        )
        location["end"] = _point_to_cds_coordinate(
            location["end"], selector_model, crossmap
        )
    else:
        location = _point_to_cds_coordinate(location, selector_model, crossmap)
    if new_variant.get("inserted"):
        new_variant["inserted"] = [
            merge_inserted_to_string(new_variant["inserted"], sequences)
        ]
    new_variant["location"] = location
    return new_variant


def reverse_start_end(variants):
    for variant in variants:
        if variant.get("location") and variant["location"]["type"] == "range":
            location = variant["location"]
            location["start"], location["end"] = location["end"], location["start"]
            location["start"]["position"] -= 1
            location["end"]["position"] -= 1


def _get_cds_into_exons(exons, cds):
    l_index = bisect.bisect_right(exons, cds[0])
    r_index = bisect.bisect_left(exons, cds[1])
    return [cds[0]] + exons[l_index:r_index] + [cds[1]]


def _location_in_same_intron(location, exons):
    start_i = bisect.bisect_right(exons, get_start(location))
    end_i = bisect.bisect_left(exons, get_end(location))
    if start_i == end_i and start_i % 2 == 0:
        return True
    else:
        return False


def _splice_site_removal(location, exons):
    start_i = bisect.bisect_right(exons, get_start(location))
    end_i = bisect.bisect_left(exons, get_end(location))
    if end_i - start_i == 1:
        return True


def _get_exons_and_cds(selector_model):
    exons = [e for l in selector_model["exon"] for e in l]
    cds = [selector_model["cds"][0][0], selector_model["cds"][0][1]]
    if selector_model.get("inverted"):
        cds[0] = exons[0]
    else:
        cds[1] = exons[-1]
    return exons, cds


def _get_exons_and_cds_2(s_m):
    exons = [e for l in s_m["exon"] for e in l]
    cds = [s_m["cds"][0][0], s_m["cds"][0][1]]
    return exons, cds


def to_exon_positions(variants, exons, cds):
    exons = _get_cds_into_exons(exons, cds)
    new_variants = []
    for variant in variants:
        if (
            variant.get("type") == "deletion_insertion"
            and variant.get("location")
            and not _location_in_same_intron(variant["location"], exons)
            and not (get_start(variant) <= exons[0] and get_end(variant) <= exons[0])
        ):
            n_v = deepcopy(variant)
            exon_s = bisect.bisect(exons, get_start(n_v))
            if exon_s % 2 == 0 and exon_s < len(exons):
                n_v["location"]["start"]["position"] = exons[exon_s]

            exon_e = bisect.bisect(exons, get_end(n_v))
            if exon_e % 2 == 0 and exon_e < len(exons):
                n_v["location"]["end"]["position"] = exons[exon_e]

            new_variants.append(n_v)

    return new_variants


def _get_splice_site_hits(variants, exons, cds):
    hits = []
    for i, variant in enumerate(variants):
        if (
            variant.get("type") == "deletion_insertion"
            and variant.get("location")
            and _splice_site_removal(
                variant["location"], _get_cds_into_exons(exons, cds)
            )
        ):
            hits.append(i)
    return hits


def reverse_variants(variants, sequences):
    reversed_variants = deepcopy(variants)
    reverse_strand_shift(reversed_variants, sequences["reference"])
    reverse_start_end(reversed_variants)
    return reversed_variants


def to_rna_protein_coordinates(variants, sequences, selector_model):
    """
    Converts the locations to cds equivalent.

    :param variants: Variants with locations in the coordinate system.
    :param sequences: Sequences with their ids as keys.
    :param selector_model: Selector model according to which
                           the conversion is performed.
    """
    exons, cds = _get_exons_and_cds(selector_model)
    crossmap = Coding(selector_model["exon"], cds, selector_model["inverted"])

    if selector_model.get("inverted"):
        variants = reverse_variants(variants, sequences)

    splice_site_hits = _get_splice_site_hits(variants, exons, cds)

    coordinate_variants = to_exon_positions(variants, exons, cds)

    cds_variants = []
    for variant in coordinate_variants:
        cds_variants.append(
            variant_to_cds_coordinate(variant, sequences, selector_model, crossmap)
        )
    return cds_variants, splice_site_hits
