import json


dict_corrections = {}
with open("cerebellum/config_data/old_regions.json") as f:
    old_regions_layer23 = json.load(f)["old_regions"]
print(old_regions_layer23)
for reg in old_regions_layer23:
    dict_corrections[reg] = [reg + 20000, reg + 30000]
# Change of id when L2 and L2/3 existed
dict_corrections[195] = [20304]
dict_corrections[747] = [20556]
dict_corrections[524] = [20582]
dict_corrections[606] = [20430]

inv_corrections = {}
for k, v in dict_corrections.items():
    for conv in v:
        inv_corrections[conv] = k
