import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file into a pandas DataFrame
df = pd.read_csv('PupilBioTest_PMP_revA.csv')

# Check the first few rows of the data
print(df.head())


from scipy.stats import ttest_ind

# Example specificity values for PMPs and CpG sites
pmp_specificity = [0.95, 0.93, 0.92, 0.96, 0.97]  # Specificity for top 5 PMPs
cpg_specificity = [0.80, 0.85, 0.79, 0.76, 0.83]  # Specificity for CpG sites

# Perform t-test
t_stat, p_value = ttest_ind(pmp_specificity, cpg_specificity)

print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")


# Combine specificity data for plotting
specificity_data = {'PMPs': pmp_specificity, 'CpG Sites': cpg_specificity}

# Create a DataFrame for easier plotting
import pandas as pd
specificity_df = pd.DataFrame(specificity_data)

# Create boxplot to compare specificity
sns.boxplot(data=specificity_df)
plt.title('Specificity Comparison: PMPs vs CpG Sites')
plt.ylabel('Specificity')
plt.show()
