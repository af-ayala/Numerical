SetFactory("OpenCASCADE");

// Parámetros
L = 1.0;      // Longitud de la viga (m)
H = 0.2;      // Altura (m)
theta = 60*Pi/180; // Ángulo en radianes
a = H / Sin(theta); // Longitud para alcanzar altura total
d = 0.1;      // Desplazamiento desde los extremos (10 cm)

// Puntos de la viga
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, H, 0, 1.0};
Point(4) = {0, H, 0, 1.0};

// Offset de grieta
x_offset = a * Cos(theta);
y_offset = a * Sin(theta);

// Puntos de grietas desplazadas desde la base
Point(5) = {d, 0, 0, 1.0};                 // Inicio grieta izquierda
Point(6) = {d + x_offset, y_offset, 0, 1.0}; // Punta grieta izquierda

Point(7) = {L - d, 0, 0, 1.0};              // Inicio grieta derecha
Point(8) = {L - d - x_offset, y_offset, 0, 1.0}; // Punta grieta derecha

// Líneas externas
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Líneas de grieta
Line(5) = {5, 6};
Line(6) = {7, 8};

// Superficie
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Malla adaptativa a lo largo de las grietas
Field[1] = Distance;
Field[1].EdgesList = {5, 6};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.002;
Field[2].SizeMax = 0.02;
Field[2].DistMin = 0.005;
Field[2].DistMax = 0.05;

Background Field = 2;
