import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('printSpline.csv', sep =';')
plt.plot(data['x'], data['y'], label='Кубический сплайн', color='blue', linewidth=2)
xx = np.linspace(data['x'].min(), data['x'].max(), 400)
plt.plot(xx, xx**3, color='orange', linewidth=2, linestyle='--', label=r'$y = x^3$')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Интерполяция кубическим сплайном')
plt.legend()
plt.grid(True)
plt.show()