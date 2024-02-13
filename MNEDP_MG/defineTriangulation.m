function [geom] = defineTriangulation(k, dVertices, ...
    dBoundary, bcBoundary, bcVertex, bcValue, checkArea, checkAngle, areaValue, ...
    angleValue, drawPlot)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%----------------------------------------------------------------------------
%
% Triangolazione di un dominio quadrata
% con condizioni di Dirichlet sul bordo
%
%----------------------------------------------------------------------------
%
%  Autore: Stefano Berrone
%  Politecnico di Torino
%
%----------------------------------------------------------------------------

%clc
%clear all

% -------------------------------
% Inserimento dei vertici
% -------------------------------

Domain.InputVertex = dVertices;


% ---------------------------------------------
% Definizione del dominio a partire dai Vertici
% ---------------------------------------------

% Dichiaro le variabili per delimitare il dominio
Domain.Boundary.Values = dBoundary;
% lato di bordo 1 dal nodo 1 al nodo 2
% lato di bordo 2 dal nodo 2 al nodo 3
% lato di bordo 3 dal nodo 3 al nodo 4
% lato di bordo 4 dal nodo 4 al nodo 1

Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio

% --------------------------------------------------
% Definizione delle condizioni al contorno a partire
% dai Vertici e dai lati di bordo
% --------------------------------------------------

% valori numerici per le condizioni al contorno
BC.Values = bcValue;

% marker delle condizioni al contorno sui bordi del dominio
% dispari -> Dirichlet; pari -> Neumann
BC.Boundary.Values = bcBoundary;
% marker dei Vertici iniziali
BC.InputVertexValues = bcVertex;
% Questi indici posso essere anche indici ai valori numerici
% contenuti nel vettore BC.Values

BC.Holes.Hole = [];
BC.Segments.Segment = [];



% --------------------------------------------
% Inserimento dei parametri di triangolazione
% --------------------------------------------

RefiningOptions.CheckArea  = checkArea;
RefiningOptions.CheckAngle = checkAngle;
RefiningOptions.AreaValue  = areaValue;
RefiningOptions.AngleValue = angleValue;
RefiningOptions.Subregions = [];


% --------------------------------------------
% Creazione della triangolazione e plottaggio
% --------------------------------------------

[geom] = mybbtr30(Domain,BC,RefiningOptions);
if drawPlot == true
    draw_grid(geom,1);
end

% --------------------------------------------------
% --------------------------------------------------


% --------------------------------------------------
% Rielaborazione dei prodotti del triangolatore
% per un piu` agevole trattamento delle condizioni
% al contorno
% --------------------------------------------------

geom.obj.P = geom.obj.P(1:geom.Nobj.N_node,:);
geom.obj.T = geom.obj.T(1:geom.Nobj.N_ele,:);
geom.obj.E = geom.obj.E(1:geom.Nobj.N_edge,:);
geom.obj.Neigh = geom.obj.Neigh(1:geom.Nobj.N_ele,:);

% --------------------------------------------------

j  = 1;
Dj = 1;
for i=1:size(geom.piv.nlist)
     if geom.piv.nlist(i)==0
        geom.piv.piv(i)=j;
        j = j+1;
     else
        geom.piv.piv(i)=-Dj;
        Dj = Dj + 1;
     end
end

% --------------------------------------------------

geom.piv.piv = transpose(geom.piv.piv);

% --------------------------------------------------

% geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
% di Dirichlet e il corrispondente marker

[X,I] = sort(geom.piv.Di(:,1));
geom.piv.Di = geom.piv.Di(I,:);

clear X I;

switch k
    case 1
        disp("Obtained P1 compatible traingulation")
    case 2
        geom = transformTriangulationP2(geom);
        disp("Obtained P2 compatible traingulation")
    otherwise
        error("Finite element method based on P" + int2str(k) + " not implemented yet.")
end
end