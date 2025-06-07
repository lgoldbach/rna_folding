

with open("phenotype_props.csv", "r") as f:
    for line in f:
        score = line.split(" ")[1]
        pre, post = score.split(".")
        comb = ".".join([pre, post[:2]])
        print(comb)