subgroup_in_order = ["uid", "max_length", "fa_threshold", "turning_angle", "tip_iteration", "smoothing", "seed_count",
                     "step_size", "min_length", "method", "otsu_threshold", "seed_plan"]
group_in_order = ["gene_condition", "sex", "strain", "age"]

# plus 1 because matlab is 1-indexed
subgroup_decoder = {"subgroup{:02d}".format(x+1): subgroup_in_order[x] for x in range(len(subgroup_in_order))}
group_decoder = {"group{}".format(x+1): group_in_order[x] for x in range(len(group_in_order))}


print(subgroup_decoder)
print(group_decoder)

MDS_column_decoder = {'group1': 'gene_condition', 'group2': 'sex', 'group3': 'strain', 'group4': 'age', 'subgroup01': 'uid', 'subgroup02': 'max_length', 'subgroup03': 'fa_threshold', 'subgroup04': 'turning_angle', 'subgroup05': 'tip_iteration', 'subgroup06': 'smoothing', 'subgroup07': 'seed_count', 'subgroup08': 'step_size', 'subgroup09': 'min_length', 'subgroup10': 'method', 'subgroup11': 'otsu_threshold', 'subgroup12': 'seed_plan'}

