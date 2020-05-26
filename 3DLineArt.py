from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

datasets = [{"x": [i, i*2/5, i*3/4], "y":[i**2, i//4, i*3],
             "z":[i, i, i], "colour": (abs(i)**(1/2)/(abs(i)+.1), abs(i)**(1/2)/(abs(i)+.1), abs(i)**(1/2)/(abs(i)+.1))} for i in range(-20, 20, 1)]

for dataset in datasets:
    ax.plot(dataset["x"], dataset["y"], dataset["z"], color=dataset["colour"])
    ax.view_init(30, 25)
plt.show()
