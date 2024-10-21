import os
import numpy as np
from scipy import io

# This python script takes array outputs of matlab script (.m) and converts them into numpy arrays (.npy)
# This is done to use the output arrays in our Convolutional Neural Network (CNN) model.

# Path to the folder containing .mat files
folder_path = r"C:\Users\erkme\OneDrive\Masaüstü\bitirme_proje_kodlar\data\mat_data\Data"

# Path to create folder which will hold the necessary numpy arrays
target_folder = r"C:\Users\erkme\OneDrive\Masaüstü\bitirme_proje_kodlar\data\Raw_data\Data\\60dB_TgnF"

Check_Data = False # Assign True if you want to observe the converted data shape.

# Function to clear a directory
def clear_directory(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)  # Remove the file
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

# Get a list of all .mat files in the folder
mat_files = [f for f in os.listdir(folder_path) if f.endswith('.mat')]

# Identify all class labels and clear corresponding folders at the start
class_labels = set()
for mat_file in mat_files:
    parts = mat_file.split('_')
    if len(parts) >= 4:
        class_labels.add(parts[3])

for class_label in class_labels:
    class_folder = os.path.join(target_folder, f"Wave_{class_label}")
    if os.path.exists(class_folder):
        clear_directory(class_folder)
    else:
        os.makedirs(class_folder, exist_ok=True)

# Loop through each .mat file
for mat_file in mat_files:
    # Extract class and sample information from the file name
    base_name, ext = os.path.splitext(mat_file)
    parts = base_name.split('_')
    if len(parts) < 4:
        print(f"Skipping file '{mat_file}' due to unexpected format")
        continue
    
    class_label = parts[2]  # Extract class label (e.g., Class_9 -> 9)
    sample_label = parts[4]  # Extract sample label (e.g., Sample_9 -> 9)
    
    # Construct the subfolder path based on the class label
    class_folder = os.path.join(target_folder, f"Wave_{class_label}")
    
    # Ensure the class folder exists
    if not os.path.exists(class_folder):
        os.makedirs(class_folder)
    
    # Load .mat file as a dictionary
    data = io.loadmat(os.path.join(folder_path, mat_file))
    
    # Convert each variable in the dictionary to a numpy array
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            # Construct the file name
            var_name = f"{base_name}_rxMatrixRawIQ"
            
            # Debug print statements to verify paths and variables
            print(f"Saving {var_name}.npy to {class_folder}")

            # Save numpy array with the constructed variable name
            np.save(os.path.join(class_folder, f"{var_name}.npy"), value)

print("Conversion completed successfully!")


# Path to the .npy file, select the path where the array you want to observe resides
npy_file_path = r"C:\Users\erkme\OneDrive\Masaüstü\bitirme_proje_kodlar\data\Raw_data\Data\\10dB\rxMatrixArrayIQ_Class_3_Sample_15_rxMatrixRawIQ.npy"


if(Check_Data == True):

    # Load the .npy file
    loaded_array = np.load(npy_file_path)

    # Display the loaded array
    print("Loaded array shape:")
    print(loaded_array.shape)