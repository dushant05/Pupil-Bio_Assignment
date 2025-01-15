import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV file into a pandas DataFrame
df = pd.read_csv('PupilBioTest_PMP_revA.csv')

# Check the first few rows of the data
print(df.head())


# Calculate the Variant Read Fraction (VRF)
df['VRF'] = df['Variant Read Count'] / df['Total Read Count']

# Display the DataFrame with VRF
print(df)

# Calculate the mean VRF for each PMP across both tissues
mean_vrf_per_pmp = df.groupby('PMP')['VRF'].mean().reset_index()

# Display the mean VRF for each PMP
print(mean_vrf_per_pmp)


# Plot the mean VRF for each PMP
plt.figure(figsize=(10, 6))
sns.barplot(x='PMP', y='VRF', data=mean_vrf_per_pmp, palette='viridis')
plt.title('Mean Variant Read Fraction (VRF) for Each PMP', fontsize=16)
plt.xlabel('PMP', fontsize=12)
plt.ylabel('Mean VRF', fontsize=12)
plt.xticks(rotation=45)
plt.show()
