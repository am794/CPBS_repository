import matplotlib.pyplot as plt

genescores_diction = []
with open('gene_scores.txt','r') as gene_scores:
    for line in gene_scores:
        columns=line.strip().split('\t')
        genescores_diction.append(columns)


loci_data={}
for row in genescores_diction:
    loci_num=row[0]
    score=row[2]

    if score != 'NA':
        if loci_num not in loci_data:
            loci_data[loci_num] = [float(score)]
        else:
            loci_data[loci_num].append(float(score))

for loci_num, scores in loci_data.items():
    mean_value = sum(scores) / len(scores)
    min_value = min(scores)
    max_value = max(scores)

    print(f"Loci {loci_num} - Mean: {mean_value}, Min: {min_value}, Max: {max_value}")

# Create a list of scores for boxplot
boxplot_data = [scores for loci_data, scores in loci_data.items()]
# Generate boxplot
plt.boxplot(boxplot_data, labels=loci_data.keys())
plt.xticks(rotation=90)
plt.xlabel('Loci Number')
plt.ylabel('Scores')
plt.title('Scores per Loci')
plt.show()