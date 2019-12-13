from pyclust import KMedoids
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

'''构造示例数据集（加入少量脏数据）'''
data1 = np.random.normal(0, 0.9, (1000, 10))
data2 = np.random.normal(1, 0.9, (1000, 10))
data3 = np.random.normal(2, 0.9, (1000, 10))
data4 = np.random.normal(3, 0.9, (1000, 10))
data5 = np.random.normal(50, 0.9, (50, 10))

data = np.concatenate((data1, data2, data3, data4, data5))

'''准备可视化需要的降维数据'''
data_TSNE = TSNE(learning_rate=100).fit_transform(data)

'''对不同的k进行试探性K-medoids聚类并可视化'''
plt.figure(figsize=(12, 8))
for i in range(2, 6):
    k = KMedoids(n_clusters=i, distance='euclidean', max_iter=1000).fit_predict(data)
    colors = ([['red', 'blue', 'black', 'yellow', 'green'][i] for i in k])
    plt.subplot(219+i)
    plt.scatter(data_TSNE[:, 0], data_TSNE[:, 1], c=colors, s=10)
    plt.title('K-medoids Resul of '.format(str(i)))
plt.show()