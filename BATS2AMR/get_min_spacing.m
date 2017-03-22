function delta = get_min_spacing(X)

X = sort(X);
dX = diff(X);

dX(dX<1e-10) = 1e99;

delta = min(abs(dX));

