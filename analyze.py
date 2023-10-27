from matplotlib import pyplot as plt
import pandas as pd
from matplotlib import gridspec
import json
# Read the data from the csv file
raw = pd.read_excel('data/results.xlsx', sheet_name='raw')
rle = pd.read_excel('data/results.xlsx', sheet_name='rle')
rle_pairing = pd.read_excel('data/results.xlsx', sheet_name='rle_pairing')
entropies = json.load(open('data/entropies.json', 'r'))
# Extract unique categories
categories = rle['Category'].unique()

# Create a separate plot for each category
for category in categories:
    plt.figure(figsize=(12, 8))  # Adjust the figure size as needed

    # Filter the DataFrame for the current category
    category_df = rle[rle['Category'] == category]
    category_df1 = rle_pairing[rle_pairing['Category'] == category]
    # Get unique error tolerances
    error_tolerances = category_df['Error Tolerance'].unique()
    plts = []
    # Create a gridspec with 2 rows and 2 columns
    gs = gridspec.GridSpec(2, 4)

    for i, error_tolerance in enumerate(error_tolerances):
        # Determine the subplot position based on the index
        if i == 0:
            subplot = plt.subplot(gs[0, :2])
        elif i == 1:
            subplot = plt.subplot(gs[0, 2:])
        elif i == 2:
            subplot = plt.subplot(gs[1, :2])
            #plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.85)
        
        plt.sca(subplot)
        plt.title(f'Error Tolerance {error_tolerance}')
        #plt.subplot(1, len(error_tolerances), i + 1)
        #plt.title(f'Error Tolerance {error_tolerance}')
        
        # Filter the DataFrame for the current error tolerance
        sub_df = category_df[category_df['Error Tolerance'] == error_tolerance]
        sub_df2 = category_df1[category_df1['Error Tolerance'] == error_tolerance]
        # Plot window size vs. interpolation time
        plt.plot(sub_df['Window Size'], sub_df['Interpolation Time'], marker='o', label='RLE', )
        plt.plot(sub_df2['Window Size'], sub_df2['Interpolation Time'], marker='o', label='RLE P')
        plts.append((sub_df['Window Size'], sub_df['Average Error'], error_tolerance))#, linestyle='--', marker='o', label='Avg Err'))
        #plt.text(.5, .05, str(f"Entropy: {round(entropies[category], 3)}"), ha='center')
        plt.legend()
        plt.plot()
        plt.xlabel('Window Size')
        plt.ylabel('Interpolation Time')

    subplot = plt.subplot(gs[1, :2])
    plt.sca(subplot)
    for p in plts:
        plt.plot(p[0], p[1], linestyle='--', marker='o', label=f'Error Tolerance {p[2]}')
    plt.legend()
    plt.plot()
    plt.xlabel('Window Size')
    plt.ylabel('Average Error')

    
    plt.suptitle(f'Data series with entropy: {round(entropies[category], 3)}')
    plt.tight_layout()
    plt.savefig(f'data/{category}.png')

plt.figure(figsize=(12, 6))
plt.bar(raw['Category'], raw['Interpolation Time'])
plt.xlabel('Category')
plt.ylabel('Interpolation Time')
plt.title('Interpolation Time on uncompressed data')
plt.savefig('data/raw.png')