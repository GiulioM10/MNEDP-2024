% Sample 11

  Domain.InputVertex = [ 0 0
    5 0
    5 5
    0 5
    1 1
    1 2
    2 1
    3 4
    4 3
    4 4];
  Domain.Boundary.Values = [1 2 3 4];
  Domain.Segments.Segment = [];
  Domain.Holes.Hole(1).Values = [5 6 7];
  Domain.Holes.Hole(2).Values = [8 9 10];

BC.Values =  [1 2 3 4 6 7];
BC.Boundary.Values = [2 1 4 3]; 
BC.Holes.Hole(1).Values = [5 5 5];
BC.Holes.Hole(2).Values = [6 6 6];
BC.Segments.Segment = [];
BC.InputVertexValues = [3 0 0 3 5 5 5 0 0 0];

RefiningOptions.CheckArea = 'Y';
RefiningOptions.CheckAngle = 'N';
RefiningOptions.AreaValue = 0.1;   
RefiningOptions.AngleValue = 30;   
RefiningOptions.Subregions = [];

[geom] = btr30(Domain,BC,RefiningOptions);
BW_NL_draw_grid (geom,1);