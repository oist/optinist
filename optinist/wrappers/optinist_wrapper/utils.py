from sklearn.preprocessing import StandardScaler


def standard_norm(X, mean, std):
    sc = StandardScaler(with_mean=mean, with_std=std)
    tX = sc.fit_transform(X)
    return tX
