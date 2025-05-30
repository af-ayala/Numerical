// Tanque de agua 3D en GMSH
// ----------------------------------------
// - Cilindro (tanque superior)
// - 4 postes inclinados (líneas)
// - Cruzamientos diagonales (opcionales)
// - Base cuadrada

// Definición de parámetros
altura_tanque = 2.0;    // Altura del cilindro
radio_tanque = 1.0;     // Radio del cilindro
altura_postes = 3.0;    // Altura de los postes
lado_base = 4.0;        // Lado de la base cuadrada
inclinacion_postes = 0.5; // Inclinación de los postes (ajustable)

// ===== 1. Geometría del tanque (cilindro) =====
// Centro del cilindro (en Z = altura_postes)
Point(1) = {0, 0, altura_postes, 1.0};
// Puntos para la base del cilindro
Point(2) = {radio_tanque, 0, altura_postes, 1.0};
Point(3) = {0, radio_tanque, altura_postes, 1.0};
Point(4) = {-radio_tanque, 0, altura_postes, 1.0};
Point(5) = {0, -radio_tanque, altura_postes, 1.0};
// Arcos del círculo base
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
// Superficie circular base
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
// Extrusión del cilindro hacia arriba
Extrude {0, 0, altura_tanque} { Surface{1}; }

// ===== 2. Postes inclinados (líneas) =====
// Puntos de la base (esquinas del cuadrado)
Point(6) = {-lado_base/2, -lado_base/2, 0, 1.0};  // Poste 1 (base)
Point(7) = {lado_base/2, -lado_base/2, 0, 1.0};   // Poste 2 (base)
Point(8) = {lado_base/2, lado_base/2, 0, 1.0};    // Poste 3 (base)
Point(9) = {-lado_base/2, lado_base/2, 0, 1.0};   // Poste 4 (base)
// Puntos superiores (con inclinación)
Point(10) = {-inclinacion_postes, -inclinacion_postes, altura_postes, 1.0};
Point(11) = {inclinacion_postes, -inclinacion_postes, altura_postes, 1.0};
Point(12) = {inclinacion_postes, inclinacion_postes, altura_postes, 1.0};
Point(13) = {-inclinacion_postes, inclinacion_postes, altura_postes, 1.0};
// Líneas de los postes
Line(5) = {6, 10};  // Poste 1
Line(6) = {7, 11};  // Poste 2
Line(7) = {8, 12};  // Poste 3
Line(8) = {9, 13};  // Poste 4

// ===== 3. Cruzamientos diagonales (opcionales) =====
// Descomentar si se quieren añadir refuerzos diagonales
/*
Line(9) = {10, 12};  // Diagonal 1 (arriba)
Line(10) = {11, 13}; // Diagonal 2 (arriba)
Line(11) = {6, 8};   // Diagonal 3 (base)
Line(12) = {7, 9};   // Diagonal 4 (base)
*/

// ===== 4. Base cuadrada (superficie) =====
Line(13) = {6, 7};  // Lado 1 (base)
Line(14) = {7, 8};  // Lado 2 (base)
Line(15) = {8, 9};  // Lado 3 (base)
Line(16) = {9, 6};  // Lado 4 (base)
Curve Loop(2) = {13, 14, 15, 16};
Plane Surface(2) = {2};

// ===== 5. Visualización y malla =====
// Estilo de visualización
Geometry.LineWidth = 2;  // Grosor de líneas
Geometry.PointSize = 3;  // Tamaño de puntos
// Tamaño de malla (ajustable)
Mesh.CharacteristicLengthMin = 0.3;
Mesh.CharacteristicLengthMax = 0.5;

// Guardar malla (opcional)
// Save "tanque_agua.msh";