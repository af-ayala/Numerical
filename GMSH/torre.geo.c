// watertower.geo
SetFactory("OpenCASCADE");

// ── PARAMETERS ─────────────────────────────────
r_cyl  = 3;    // radius of the stem (cylinder)
h_cyl  = 30;   // height of the stem

r_bot  = r_cyl;// bottom radius of the cone (matches top of stem)
r_top  = 12;   // top radius of the cone
h_cone = 10;   // height of the cone

// ── GEOMETRY ───────────────────────────────────
// 1) Cylinder (stem)
Cylinder(1) = {
  0, 0, 0,    // base center
  0, 0, h_cyl,// axis vector (height along Z)
  r_cyl       // radius
};

// 2) Truncated cone (tank) sitting on top of cylinder
//    Cone(tag) = {XO, YO, ZO, DX, DY, DZ, R1, R2}
Cone(2) = {
  0, 0, h_cyl,  // base center sits at top of cylinder
  0, 0, h_cone, // axis vector of length = h_cone
  r_bot,        // base radius
  r_top         // top radius
};

// 3) Fuse into one single volume
BooleanUnion{
  Volume{1};  // cylinder
  Volume{2};  // cone
  Delete;     // delete the originals
}{
  Volume{3};  // union result
}

// ── MESH SETTINGS ──────────────────────────────
Mesh.CharacteristicLengthMin = 1.0;
Mesh.CharacteristicLengthMax = 2.0;
