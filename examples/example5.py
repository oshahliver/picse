"""This is an example of using pre-generated data sets to create new model calibrations
using various approaches including simple regression models and neural networks
"""

from picse.utils.file_tools import internal_data as intdat
import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from sklearn.model_selection import train_test_split
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler

filepath = "/home/os18o068/Documents/PHD/Projects/pics_external_data/training_sets/training_1.csv"
df, meta = intdat.read_sample(filepath)
# df = df.sample(n = 200000)

# extract the features and labels from the data frame
X = df[
    [
        "mass",
        "pres_surface",
        "si_number_mantle",
        "mg_number",
        "fe_number_mantle",
        "ocean_mass_fraction",
        "temp_surface",
    ]
]
y = df[["pres_center"]]

# take the log of mass and pres_center
X["mass"] = np.log(X["mass"])
X["pres_surface"] = np.log(X["pres_surface"])
X["temp_surface"] = np.log(X["temp_surface"])
# y["temp_center"] = np.log(y["temp_center"])
y["pres_center"] = np.log(y["pres_center"])

# split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# normalize the features and target variables
# to mean = 0 and std = 1
scaler = StandardScaler()
X_train_norm = scaler.fit_transform(X_train)
X_test_norm = scaler.fit_transform(X_test)

# set up parameter ranges
num_layers = [6, 7, 8]
layer_size = [128, 256, 512]
activation = ["relu"]

# iterate over all combinations of parameters
best_loss = float("inf")
pbar1 = tqdm(total=len(num_layers) * len(layer_size) * len(activation))
for nl in num_layers:
    for ls in layer_size:
        for act in activation:
            # define the neural network model
            model = Sequential()
            model.add(Dense(ls, input_shape=(7,), activation=act))
            for i in range(nl - 1):
                model.add(Dense(ls, activation=act))
            model.add(Dense(2, activation="linear"))

            # compile the model
            model.compile(optimizer="adam", loss="mape")

            # train the model
            n_iter = 50
            pbar2 = tqdm(total=n_iter, leave=False, desc=f"NL={nl}, LS={ls}, ACT={act}")
            for i in range(n_iter):
                model.fit(X_train_norm, y_train, epochs=1, batch_size=32, verbose=0)
                pbar2.update(1)

            pbar1.update(1)

            # evaluate the model on the test set
            test_loss = model.evaluate(X_test_norm, y_test, verbose=0)
            print(
                f"Number of layers: {nl}, Layer size: {ls}, Activation: {act}, Test loss: {test_loss:.3f}"
            )

            # check if this is the best model so far
            if test_loss < best_loss:
                best_loss = test_loss
                best_model = model
                best_parameters = (nl, ls, act)

# use the best model to make predictions on the test set
predictions = best_model.predict(X_test_norm)

print(f"Best model test loss: {best_loss:.3f}")
print(
    f"Best parameters: Number of layers: {best_parameters[0]}, Layer size: {best_parameters[1]}, Activation: {best_parameters[2]}"
)
print(predictions)
