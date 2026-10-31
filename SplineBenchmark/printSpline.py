import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('printSpline.csv', sep =';')
plt.plot(data['x'], data['y'], label='Кубический сплайн', color='blue', linewidth=2)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Интерполяция кубическим сплайном')
plt.legend()
plt.grid(True)
plt.show()