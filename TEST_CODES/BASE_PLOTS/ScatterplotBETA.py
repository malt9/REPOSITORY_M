import kagglehub
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Download latest version
path = kagglehub.dataset_download("denkuznetz/housing-prices-regression")
dataset_path = path+"/real_estate_dataset.csv" 
dataset = pd.read_csv(dataset_path)
# dataset.info()

#Prepare useful data
cols = list(dataset.columns)
cols.remove("ID")
cols.remove("Price")



fig, axs = plt.subplots(nrows=2, ncols=5)

for ax, col in zip(axs.flat, cols):
    modelo = linregress(dataset[col], dataset["Price"])
    ax.scatter(dataset[col], dataset["Price"])
    print(f"Coeficiente de Correlación [{col}]:",modelo.rvalue)
    print("intercepto:",modelo.intercept)
    print("pendiente:",modelo.slope)

plt.show()