import copy

from normalizer.util import (
    roll,
    get_start,
    get_end,
    get_location_length,
    get_inserted_length,
)


def is_deletion(delins_variant):
    if (delins_variant.get("inserted") is None) or (
        len(delins_variant.get("inserted")) == 0
    ):
        return True
    for inserted in delins_variant["inserted"]:
        if get_location_length(inserted["location"]) > 0:
            return False
    return True

def update_inserted_with_sequences(inserted, sequences):
    for insert in inserted:
        if insert["source"] == "observed":
            insert["sequence"] = sequences["observed"][
                get_start(insert["location"]) : get_end(insert["location"])
            ]


def seq2repeats(input):
    for i in range(1, len(input) + 1):
        for j in range(1, len(input) + 1):
            if i * j == len(input):
                if input[:j] * i == input:
                    yield(input[:j], i)
                break

def seq_present_before(observed, ins_seq, start, end):

    for repeat, i in seq2repeats(ins_seq):
        if observed[start - len(repeat): end] == repeat:
            print(f"found repeat: {repeat}")
            return repeat, i + 1

    return "", 0

def de_to_hgvs(variants, sequences=None):
    """
    Convert the variants to an HGVS format (e.g., a deletion insertion
    of one nucleotide is converted to a substitution).
    """
    new_variants = []
    o_index = 0
    if len(variants) == 1 and variants[0].get("type") == "equal":
        return [copy.deepcopy(variants[0])]
    for variant in variants:
        if variant.get("type") == "equal":
            o_index += get_location_length(variant["location"])
        elif variant.get("type") == "inversion":
            o_index += get_location_length(variant["location"])
            new_variants.append(copy.deepcopy(variant))
        elif variant.get("type") == "deletion_insertion":
            if get_start(variant["location"]) == get_end(variant["location"]):
                # delins_to_insertion
                shift5, shift3 = roll(
                    sequences["observed"],
                    o_index + 1,
                    o_index + get_inserted_length(variant["inserted"]),
                )
                ins_seq = sequences["observed"][
                    o_index : o_index + get_inserted_length(variant["inserted"])
                ]
                o_index += shift3
                new_variant = copy.deepcopy(variant)
                ins_length = get_location_length(variant["inserted"][0]["location"])
                new_variant["location"]["start"]["position"] += shift3
                new_variant["location"]["end"]["position"] += shift3

                repeat_seq, repeat_number = seq_present_before(sequences["observed"], ins_seq,
                                   get_start(variant["location"]),
                                   get_end(variant["location"]))

                if (
                    sequences["observed"][
                        get_start(variant["location"])
                        - ins_length : get_end(variant["location"])
                    ]
                    == ins_seq
                ):
                    new_variant["type"] = "duplication"
                    new_variant["location"]["start"]["position"] = (
                        get_start(new_variant["location"]) - ins_length
                    )
                elif len(repeat_seq) > 0:
                    new_variant["type"] = "repeat"
                    new_variant["inserted"][0]["repeat_number"] = {"value": repeat_number}
                    new_variant["inserted"][0]["location"]["end"]["position"] = \
                        get_start(new_variant["location"]) + len(repeat_seq)
                    variant["inserted"][0]["source"] = "description"
                    new_variant["location"]["start"]["position"] = (
                        get_start(new_variant["location"]) - len(repeat_seq)
                    )
                else:
                    new_variant["type"] = "insertion"

                update_inserted_with_sequences(new_variant["inserted"], sequences)
                new_variants.append(new_variant)

            elif is_deletion(variant):
                # delins_to_deletion
                new_variant = {
                    "type": "deletion",
                    "source": "reference",
                    "location": copy.deepcopy(variant["location"]),
                }
                shift5, shift3 = roll(
                    sequences["reference"],
                    new_variant["location"]["start"]["position"] + 1,
                    new_variant["location"]["end"]["position"],
                )
                new_variant["location"]["start"]["position"] += shift3
                new_variant["location"]["end"]["position"] += shift3
                new_variants.append(new_variant)

            elif len(variant["inserted"]) == 1:
                if (
                    variant["inserted"][0]["location"] == variant["location"]
                    and variant["inserted"][0]["source"] == "reference"
                ):
                    if variant["inserted"][0].get("inverted") is True:
                        # delins_to_inversion
                        new_variant = {
                            "type": "inversion",
                            "source": "reference",
                            "location": copy.deepcopy(variant["location"]),
                        }
                        new_variants.append(new_variant)
                elif (
                    get_location_length(variant["location"])
                    == get_location_length(variant["inserted"][0]["location"])
                    == 1
                ):
                    # delins_to_substitution
                    new_variant = copy.deepcopy(variant)
                    new_variant["type"] = "substitution"
                    update_inserted_with_sequences(new_variant["inserted"], sequences)
                    new_variant["deleted"] = [{
                        "sequence": sequences["reference"][
                            get_start(new_variant["location"]) : get_end(
                                new_variant["location"]
                            )
                        ],
                        "source": "reference_location",
                    }]
                    new_variants.append(new_variant)
                else:
                    # delins_to_deletion_insertion
                    new_variant = copy.deepcopy(variant)
                    update_inserted_with_sequences(new_variant["inserted"], sequences)
                    new_variants.append(new_variant)
            else:
                # TODO: Check if this is a deletion insertion and if the
                # sequence should or not be merged.
                # TODO: Checkif a shift should be performed also here.
                new_variant = copy.deepcopy(variant)
                update_inserted_with_sequences(new_variant["inserted"],
                                               sequences)
                new_variants.append(new_variant)

            o_index += get_inserted_length(variant["inserted"])

        else:
            raise ValueError(
                "Unexpected variant type: '{}'.".format(variant.get("type"))
            )

    return new_variants
