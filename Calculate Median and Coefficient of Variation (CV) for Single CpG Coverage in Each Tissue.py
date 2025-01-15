import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV file into a pandas DataFrame
df = pd.read_csv('PupilBioTest_PMP_revA.csv')

# Check the first few rows of the data
print(df.head())


# Group by tissue and calculate the median coverage
median_coverage = df.groupby('Tissue')['Coverage'].median()
print("Median Coverage by Tissue:")
print(median_coverage)

# Group by tissue and calculate the mean and standard deviation
mean_coverage = df.groupby('Tissue')['Coverage'].mean()
std_coverage = df.groupby('Tissue')['Coverage'].std()

# Calculate the Coefficient of Variation (CV) for each tissue
cv = (std_coverage / mean_coverage) * 100
print("Coefficient of Variation (CV) by Tissue:")
print(cv)

# Combine the results into a single DataFrame
summary_df = pd.DataFrame({'Median Coverage': median_coverage, 'CV (%)': cv})

# Print the summary DataFrame
print("Summary of Median and CV:")
print(summary_df)


# Set the style for the plot
sns.set(style="whitegrid")

# Create a histogram to visualize the distribution of CpG coverage
plt.figure(figsize=(8, 6))
sns.histplot(df['Coverage'], kde=True, bins=30, color='skyblue', edgecolor='black')
plt.title('Distribution of CpG Coverage', fontsize=16)
plt.xlabel('Coverage', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.show()

# Create a boxplot for coverage by tissue type
plt.figure(figsize=(8, 6))
sns.boxplot(x='Tissue', y='Coverage', data=df)
plt.title('CpG Coverage Distribution by Tissue')
plt.xlabel('Tissue')
plt.ylabel('CpG Coverage')
plt.show()

# Create a violin plot to visualize the coverage distribution by tissue type
plt.figure(figsize=(8, 6))
sns.violinplot(x='Tissue', y='Coverage', data=df, palette='Set2')
plt.title('Violin Plot of CpG Coverage Across Tissues', fontsize=16)
plt.xlabel('Tissue', fontsize=12)
plt.ylabel('CpG Coverage', fontsize=12)
plt.show()

# Assuming another column, e.g., 'Gene_Expression', is available
# Scatter plot of coverage vs gene expression (if applicable)
plt.figure(figsize=(8, 6))
sns.scatterplot(x=df['Coverage'], y=df['Gene_Expression'], hue=df['Tissue'], palette='Set1')
plt.title('Scatter Plot of CpG Coverage vs Gene Expression', fontsize=16)
plt.xlabel('CpG Coverage', fontsize=12)
plt.ylabel('Gene Expression', fontsize=12)
plt.show()


# Create a density plot for coverage values
plt.figure(figsize=(8, 6))
sns.kdeplot(data=df, x='Coverage', hue='Tissue', fill=True, common_norm=False, palette='Set2')
plt.title('Density Plot of CpG Coverage Across Tissues', fontsize=16)
plt.xlabel('CpG Coverage', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()


# Overlay histogram and KDE plot for coverage distribution
plt.figure(figsize=(8, 6))
sns.histplot(df['Coverage'], kde=True, bins=30, color='skyblue', edgecolor='black', stat='density')
plt.title('Overlay Histogram and KDE of CpG Coverage', fontsize=16)
plt.xlabel('CpG Coverage', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()
