from confusion_matrix import df
import json

# Checks whether the new code produces the same confusion matrix as the old code using artificial data.
with open("artificial_confusion_matrix.json", "r") as file:
    correct_df = json.loads(file.read())

if df == correct_df:
    print("test passed")
else:
    print("test failed")