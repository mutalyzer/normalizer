from .util import get_start, get_end, update_position, sort_variants
from .converter import get_mol_type
from .references import get_transcripts_ids


def are_sorted(variants):
    """
    Check if the provided variants list is sorted.
    """
    current_position = 0
    for variant in variants:
        if get_start(variant['location']) < current_position:
            return False
        current_position = get_start(variant['location'])
    return True


def is_overlap(variants):
    """
    Check whether there is overlap between the variants.

    TODO: Add support for fuzzy (uncertain) locations.
    """
    sorted_variants = sort_variants(variants)

    current_end_position = 0
    for variant in sorted_variants:
        if get_start(variant['location']) < current_end_position:
            return True
        else:
            current_end_position = get_end(variant['location'])
    return False


def check_location_start_end(location, indexing='internal'):
    if indexing == 'internal':
        return get_start(location) <= get_end(location)
    elif indexing == 'hgvs':
        if location['type'] == 'point':
            return True
        else:
            return get_start(location) < get_end(location)


def check_variant_start_end(variant, reference=None, indexing='internal'):
    if indexing == 'internal':
        return get_start(variant['location']) <= get_end(variant['location'])
    elif indexing == 'hgvs':
        return get_start(variant['location']) < get_end(variant['location'])


def check_substitution(variant, reference=None):
    pass


def check_deletion(variant, reference=None):
    pass




def check_semantics(variants, reference=None, indexing='internal'):
    start_end = []
    for variant in variants:
        start_end.append(check_variant_start_end(variant))
    output = {'sorted': are_sorted(variants),
              'overlap': is_overlap(variants),
              'start_end': start_end}

    return output


def is_fuzzy_point(point_location):
    if point_location.get('uncertain'):
        return True
    if point_location.get('offset') and \
            point_location['offset'].get('uncertain'):
        return True
    return False


def is_fuzzy_range(range_location):
    if range_location.get('uncertain'):
        return True
    if is_fuzzy_point(range_location['start']):
        return True
    if is_fuzzy_point(range_location['end']):
        return True
    return False


def is_fuzzy_location(location):
    if location['type'] == 'range':
        return is_fuzzy_range(location)
    if location['type'] == 'point':
        return is_fuzzy_point(location)


def check_for_fuzzy(variants):
    for variant in variants:
        if variant.get('location') and is_fuzzy_location(variant['location']):
            return True
        if variant.get('inserted'):
            for inserted in variant['inserted']:
                if inserted.get('location') and \
                        is_fuzzy_location(inserted['location']):
                    return True
    return False


def is_intronic_point(point_location):
    if point_location.get('offset') and point_location['offset']['value'] != 0:
        return True
    return False


def is_intronic_range(range_location):
    if is_intronic_point(range_location['start']):
        return True
    if is_intronic_point(range_location['end']):
        return True


def is_intronic_location(location):
    if location['type'] == 'range':
        return is_intronic_range(location)
    elif location['type'] == 'point':
        return is_intronic_point(location)


def check_intronic_positions(variants):
    for variant in variants:
        if variant.get('location') and \
                is_intronic_location(variant['location']):
            return True
        if variant.get('inserted'):
            for inserted in variant['inserted']:
                if inserted.get('location') and \
                        is_intronic_location(inserted['location']):
                    return True
    return False


def is_out_of_range(location, sequence):
    if location['type'] == 'point':
        if not 0 <= location['position'] <= len(sequence):
            return True
    if location['type'] == 'range':
        if not 0 <= location['start']['position'] <= len(sequence):
            return True
        if not 0 <= location['end']['position'] <= len(sequence):
            return True


def check_out_of_range(variants, sequences):
    """
    Check if the provided locations are within the reference length.
    """
    for variant in variants:
        if variant.get('location') and \
                is_out_of_range(variant['location'], sequences['reference']):
            return True
        if variant.get('inserted'):
            for inserted in variant['inserted']:
                if inserted.get('location'):
                    if isinstance(inserted['source'], dict):
                        if is_out_of_range(
                                inserted['location'],
                                sequences[inserted['source']['id']]):
                            return True
                    elif inserted.get('source') == 'reference':
                        if is_out_of_range(
                                inserted['location'],
                                sequences[inserted['source']]):
                            return True
    return False


def check_description_sequences(variants, sequence):
    """
    Check if there is a match between the actual sequence and the
    user provided ones. Examples:
    - NG_029724.1:g.10_11delGT: check if GT is located at 10_11.
    """
    for variant in variants:
        if variant.get('deleted'):
            for deleted in variant['deleted']:
                if deleted.get('sequence') and \
                        sequence[variant['location']['start']['position']:
                        variant['location']['end']['position']] != \
                        deleted['sequence']:
                    return True
                if deleted.get('length') and \
                        variant['location']['end']['position'] - \
                        variant['location']['start']['position'] != \
                        deleted['length']['value']:
                    return True

        if variant['type'] == 'duplication':
            if variant.get('inserted'):
                for inserted in variant['inserted']:
                    if sequence[
                       variant['location']['start']['position']:
                       variant['location']['end']['position']] != \
                            inserted['sequence']:
                        return True
    return False


def check_coordinate_system(description_model, references):
    """
    Check if there is a match between the provided coordinate system the
    reference type.
    - A c. coordinate system can be used with a genomic reference only if a
    selector is provided.
    """
    reference = references[description_model['reference']['id']]
    if description_model.get('coordinate_system') == 'c':
        if get_mol_type(reference) == 'genomic DNA' and \
                (description_model['reference'].get('selector') is None):
            raise Exception('No selector mentioned for c. with a genomic '
                            'reference.')


def get_location_length(variant):
    return get_end(variant['location']) - get_start(variant['location'])


def check_location_length_is_one(variant):
    if get_location_length(variant) != 1:
        raise Exception('Not consecutive positions in insertion.')


def check_length_in_inserted(variant):
    for inserted in variant['inserted']:
        if inserted.get('length'):
            raise Exception('Length in inserted not supported.')


def check_positions(variant):
    if variant['type'] == 'insertion':
        check_location_length_is_one(variant)
        check_length_in_inserted(variant)
    if variant['type'] == 'substitution':
        check_length_in_inserted(variant)


def is_unsupported_variant_type(variant):
    if variant['type'] not in ['substitution', 'deletion', 'duplication',
                               'insertion', 'inversion', 'conversion',
                               'deletion_insertion']:
        return True
    else:
        return False


def validate_variants(variants, sequences):
    for variant in variants:
        if is_unsupported_variant_type(variant):
            raise Exception('Variant type not supported.')
        check_positions(variant)
