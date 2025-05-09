// Dominio cuadrado 1x1 m con malla 100x100 (h = 0.01 m)
lc = 0.01;  // Tamaño característico

// Puntos del cuadrado
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Líneas del contorno
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Superficie
Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

// Malla 2D

// Mesh 2;
// SetFactory("OpenCASCADE");