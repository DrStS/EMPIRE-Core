GiD Post Results File 1.0

GaussPoints "GP_triangles" ElemType Triangle "mesh_triangles"
  Number Of Gauss Points: 1
  Natural Coordinates: internal
end gausspoints

GaussPoints "GP_quads" ElemType Quadrilateral "mesh_quads"
  Number Of Gauss Points: 1
  Natural Coordinates: internal
end gausspoints

Result	"nodal data"	"EMPIRE_CoSimulation"	5	Vector	OnNodes
Values
	3	3	0	0
	1	1	0	0
	5	5	0	0
	2	2	0	0
	6	6	0	0
End Values

Result	"elemental data"	"EMPIRE_CoSimulation"	6	Vector	OnGaussPoints	"GP_triangles"
Values
	2	20	0	0
End Values

Result	"elemental data"	"EMPIRE_CoSimulation"	6	Vector	OnGaussPoints	"GP_quads"
Values
	1	10	0	0
End Values

