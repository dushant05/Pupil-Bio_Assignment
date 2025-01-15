import pandas as pd

# Load the CSV file into a pandas DataFrame
df = pd.read_csv('PupilBioTest_PMP_revA.csv')

# Check the first few rows of the data
print(df.head())



from sklearn.preprocessing import StandardScaler

# Assuming df contains the expression data (excluding 'PMP' column)
scaler = StandardScaler()
normalized_data = scaler.fit_transform(df.drop('PMP', axis=1))
df_normalized = pd.DataFrame(normalized_data, columns=df.columns[1:], index=df['PMP'])

from scipy import stats

# Assuming tissue_1 and tissue_2 are expression arrays for a specific PMP
t_stat, p_val = stats.ttest_ind(df['Tissue #1 Expression'], df['Tissue #2 Expression'])

from statsmodels.stats.multitest import multipletests

# Assuming p_values is a list of p-values from the t-tests
_, corrected_pvals, _, _ = multipletests(p_values, method='fdr_bh')




from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

# Create labels: 1 for Tissue #1, 0 for others
y = (df['Tissue'] == 'Tissue #1').astype(int)
X = df.drop(['PMP', 'Tissue'], axis=1)

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize the Random Forest model
rf = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
rf.fit(X_train, y_train)

# Predict on the test set
y_pred = rf.predict(X_test)

# Get the classification report
print(classification_report(y_test, y_pred))


feature_importance = rf.feature_importances_
importance_df = pd.DataFrame({
    'PMP': X.columns,
    'Importance': feature_importance
}).sort_values(by='Importance', ascending=False)

print(importance_df)


probs = rf.predict_proba(X_test)[:, 1]



threshold = 0.9  # Only select PMPs with > 90% confidence
high_confidence_pmp = probs >= threshold


confidence_scores = (corrected_pvals < 0.05) & (probs > 0.8)  # Define high-confidence PMPs


final_pmp_list = pd.DataFrame({
    'PMP': X.columns,
    'P-value': corrected_pvals,
    'Confidence Score': confidence_scores
}).sort_values(by='Confidence Score', ascending=False)

print(final_pmp_list)




from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

fpr, tpr, thresholds = roc_curve(y_test, probs)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()
