#!/bin/env python3

def pos_to_region(pos, starts, ends, names):
    '''Takes a list of positions and regions, and returns lists of region names that contain the position.'''

    result = dict()

    # regions must be storted by start
    s_starts, s_names = [list(i) for i in zip(*[(x, y) for x, y in sorted(zip(starts, names))])]
    e_ends, e_names = [list(i) for i in zip(*[(x, y) for x, y in sorted(zip(ends, names))])]
    
    if type(pos) != int:
        s_pos = sorted(pos)
    else:
        return 'Not implemented for single position'

    open_pos = []
    next_start = 0
    next_end = 0

    for p in s_pos:

        if next_start < len(s_starts):
            while p >= s_starts[next_start]:
                open_pos.append(s_names[next_start])
                next_start += 1
                if next_start >= len(s_starts): break

        if next_end < len(e_ends):
            while p > e_ends[next_end]:
                open_pos.remove(e_names[next_end])
                next_end += 1
                if next_end >= len(e_ends): break

        result[p] = list(open_pos)
        
        if next_end > len(e_ends): break
        
    return [','.join(result[p]) for p in pos] # return results in the provided order

def region_to_pos(pos, starts, ends, names):
    '''Takes a list of positions and regions, and returns lists of positions that are contained in each region.'''
    result = {n: [] for n in names}

    # regions must be storted by start
    s_starts, s_names = [list(i) for i in zip(*[(x, y) for x, y in sorted(zip(starts, names))])]
    e_ends, e_names = [list(i) for i in zip(*[(x, y) for x, y in sorted(zip(ends, names))])]
    s_pos = sorted(pos)
    
    open_pos = []
    next_start = 0
    next_end = 0

    for p in s_pos:

        if next_start < len(s_starts):
            while p >= s_starts[next_start]:
                open_pos.append(s_names[next_start])
                next_start += 1
                if next_start >= len(s_starts): break

        if next_end < len(e_ends):
            while p > e_ends[next_end]:
                open_pos.remove(e_names[next_end])
                next_end += 1
                if next_end >= len(e_ends): break
        
        for n in open_pos:
            result[n] += [str(p)]

        if next_end > len(e_ends): break
    
    return [','.join(result[n]) for n in names] # return results in the provided order

def region_to_n_mut(pos, starts, ends, names):
    '''Takes a list of positions and regions, and returns lists of positions that are contained in each region.'''
    result = {n: 0 for n in names}

    # regions must be storted by start
    s_starts, s_names = [list(i) for i in zip(*[(x, y) for x, y in sorted(zip(starts, names))])]
    e_ends, e_names = [list(i) for i in zip(*[(x, y) for x, y in sorted(zip(ends, names))])]
    s_pos = sorted(pos)
    
    open_pos = []
    next_start = 0
    next_end = 0

    for p in s_pos:

        if next_start < len(s_starts):
            while p >= s_starts[next_start]:
                open_pos.append(s_names[next_start])
                next_start += 1
                if next_start >= len(s_starts): break

        if next_end < len(e_ends):
            while p > e_ends[next_end]:
                open_pos.remove(e_names[next_end])
                next_end += 1
                if next_end >= len(e_ends): break
        
        for n in open_pos:
            result[n] += 1

        if next_end > len(e_ends): break
    
    return [result[n] for n in names] # return results in the provided order

def cut_overlapping_regions(df, region_names = 'name', starts_names = 'start', ends_names = 'end', sep = '+', hole_sep = '_to_', genome_size = 4641652, include_holes = False):
    '''Takes a dataframe with regions (as start and end) and transform it in a dataframe with the overlapping regions cut into their own section.'''
    
    import pandas as pd
    
    # Convert position list to event list
    event_lines = []

    for l in range(len(df)):

        line = df.iloc[l]
        name = line[region_names]
        start = line[starts_names]
        end = line[ends_names] + 1

        event_lines += [[name, 'start', start]]
        event_lines += [[name, 'end', end]]

    events = pd.DataFrame(event_lines, columns = ['region', 'event', 'pos'])

    # Sort events
    events.sort_values(['pos', 'event'], ascending = [True, False], inplace = True)

    # Convert events to sections
    open_genes = [] # genes that are overlapping the current cursor position
    last_gene = 'genome_start' # last gene to be closed, for intergenic regions
    last_pos = 1
    regions = []

    for l in range(len(events)):

        # Is the current region a hole?
        if not open_genes:

            # Open a new gene
            assert events.iloc[l].event == 'start'
            line = [last_gene + hole_sep + events.iloc[l].region, last_pos, events.iloc[l].pos - 1]
            open_genes += [events.iloc[l].region]
            last_gene = events.iloc[l].region

            if last_pos != events.iloc[l].pos and include_holes:
                regions += [line]

        else:
            # Are we adding a new overlap?
            if events.iloc[l].event == 'start':

                # Write the genes that were already open (not including the new gene)
                line = [sep.join(open_genes), last_pos, events.iloc[l].pos - 1]
                open_genes += [events.iloc[l].region]
                last_gene = events.iloc[l].region

            # Are we closing a gene?
            elif events.iloc[l].event == 'end':
                # Write the genes that were already open (current gene should be here)
                assert events.iloc[l].region in open_genes
                line = [sep.join(open_genes), last_pos, events.iloc[l].pos - 1]

                open_genes.remove(events.iloc[l].region)

            # Add to the list of regions, unless it's a zero-length region
            if last_pos != events.iloc[l].pos:
                regions += [line]

        last_pos = events.iloc[l].pos

    # Add a last line to close the circular chromosome
    if include_holes:
        regions += [[last_gene + '_to_genome_end', last_pos, genome_size]]

    # Aggregate into a DF
    return pd.DataFrame(regions, columns = [region_names, starts_names, ends_names])