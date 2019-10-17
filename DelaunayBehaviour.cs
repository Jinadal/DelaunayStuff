using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Troll3D.Components;
using Troll3D;
using SharpDX;

using DelaunayTriangularisation.WingedEdge;

namespace DelaunayTriangularisation
{
    public class DelaunayBehaviour : Behaviour
    {
        public MeshRenderer meshRenderer;

        public void CleanVoronoi()
        {
            for ( int i = 0; i < VoronoiPoints.Count; i++ )
            {
                Scene.CurrentScene.RemoveRenderable( ( MeshRenderer )VoronoiPoints[i].GetComponent( ComponentType.MeshRenderer ) );
            }
            VoronoiPoints.Clear();
            ( ( LineRenderer )VoronoiLines.GetComponent( ComponentType.LineRenderer ) ).Display = false;
        }

        public void BuildVoronoi(WingedEdgeMesh mesh)
        {
             List<StandardVertex> lines = new List<StandardVertex>();

			// We start by building the vertices corresponding to the face's center of the triangulation who are node of the Voronoi diagram
             CleanVoronoi();
            
            for ( int i = 0; i < mesh.Faces.Count; i++ )
            {
                // We get back the face's center (barycentre)

                List<VertexWE> vertices = mesh.GetFaceVertices( mesh.Faces[i] );

                Vector3 center = vertices[0].Position + vertices[1].Position + vertices[2].Position;
                center = center / 3.0f;

                Troll3D.Entity entity = new Troll3D.Entity( Entity );
                entity.transform_.SetPosition(
                    center.X,
                    center.Y,
                    0.0f );

                entity.transform_.SetScale( 0.02f, 0.02f, 1.0f );

                MaterialDX11 material = new MaterialDX11();
                material.SetMainColor( 0.0f, 1.0f, 1.0f, 1.0f );
                MeshRenderer meshrenderer = entity.AddComponent<MeshRenderer>();
                meshrenderer.material_ = material;
                meshrenderer.model_ = Quad.GetMesh();

                // We get back the neighbours

                List<FaceWE> neighbours = mesh.GetFaceNeighbours( mesh.Faces[i] );

                for ( int j = 0; j < neighbours.Count; j++ )
                {
                    // On récupère le centre de la face (barycentre)

                    List<VertexWE> verticesNeighbour = mesh.GetFaceVertices( neighbours[j] );

                    Vector3 centerNeighbour = verticesNeighbour[0].Position + verticesNeighbour[1].Position + verticesNeighbour[2].Position;
                    centerNeighbour = centerNeighbour / 3.0f;

                    lines.Add( new StandardVertex( center ) );
                    lines.Add( new StandardVertex( centerNeighbour ) );
                }

                VoronoiPoints.Add( entity );
            }

            if ( lines.Count > 0 )
            {
                ( ( LineRenderer )VoronoiLines.GetComponent( ComponentType.LineRenderer ) ).Display = true;
                ( ( LineRenderer )VoronoiLines.GetComponent( ComponentType.LineRenderer ) ).Vertices = lines;
                ( ( LineRenderer )VoronoiLines.GetComponent( ComponentType.LineRenderer ) ).UpdateRenderer();
            }
            
        }

        public List<Entity> VoronoiPoints = new List<Entity>();
        public Entity VoronoiLines;

        public override void Initialize()
        {
            VoronoiLines = new Entity();
            LineRenderer lr =  VoronoiLines.AddComponent<LineRenderer>();
            lr.material_ = new MaterialDX11("vDefault.cso","pUnlit.cso","gDefaultLine.cso");
            lr.material_.SetMainColor( 0.0f, 0.0f, 1.0f, 1.0f );

            meshRenderer = Entity.AddComponent<MeshRenderer>();
            meshRenderer.material_ = new MaterialDX11();
        }

        /// <summary>
		/// This function make sure point P is insibe triangle ABC. 
		/// We compute the matrix determinant then we check if the result is > 0 or not.
		/// ABC must be a triangle in clockwise order
        /// </summary>
        /// <returns></returns>
        public bool IsPointInTriangle( Vector3 P, Vector3 A, Vector3 B, Vector3 C )
        {
            float [] matvalues = new float[16]
            {
                A.X, A.Y, (float)(Math.Pow(A.X,2)+Math.Pow(A.Y,2)),1,
                B.X, B.Y, (float)(Math.Pow(B.X,2)+Math.Pow(B.Y,2)),1,
                C.X, C.Y, (float)(Math.Pow(C.X,2)+Math.Pow(C.Y,2)),1,
                P.X, P.Y, (float)(Math.Pow(P.X,2)+Math.Pow(P.Y,2)),1
            };
            return ((new Matrix(matvalues)).Determinant()>0);
        }

        /// <summary>
		/// To find the center of the circumscribed circle of the triangle, we need to find
		/// the intersection point of the 3 perpendicular bisector of the triangle
        /// </summary>
        /// <returns></returns>
        public Vector3 GetCircleCenter( Vector3 A, Vector3 B, Vector3 C )
        {
            Vector3 AB = ( B - A );
            Vector3 BC = ( C - B );
            Vector3 CA = ( A - C );

            // Je récupère un produit vectoriel que je pourai ensuite utiliser pour
            // trouver les médiatrices

			// Using Cross product to find perpendicular bisectors
            Vector3 cross = Vector3.Cross(AB,BC);
            cross.Normalize();

            // Vertex from where the perpendicular bisectors comes off
            Vector3 M = A + (AB / 2.0f);
            Vector3 U = Vector3.Cross( AB, cross );
            U.Normalize();

            Vector3 N = B + (BC / 2.0f);
            Vector3 V = Vector3.Cross( BC, cross );
            V.Normalize();

            Vector3 result = new Vector3() ;
            MathStuffs.SegmentIntersect( M, M + U, N, N + V, ref result );

			// Now looking for the intersection point between these 2 line segments

            // Pour se faire, on va résoudre l'équation à 2 inconnues utilisant l'équation
            //// paramétriques des 2 droites

            //float s = M.Y + ( ( N.X * U.Y ) / U.X ) - ( ( M.X * U.Y )/ U.X  ) - N.Y;
            //float denominateur = U.Y  - (( V.X * U.Y ) / U.X);
            //s = s / denominateur;

            //float x = N.X + s * V.X;
            //float y = N.Y + s * V.Y;

            return new Vector3( result.X, result.Y, 0.0f );
        }

        /// <summary>
        /// Retourne le rayon d'un cercle circonscrit à un triangle de centre P
        /// </summary>
        /// <returns></returns>
        public float GetCircleRadius( Vector3 A, Vector3 B, Vector3 C, Vector3 P )
        {
            float maxradius = -1.0f;

            if ( ( P - A ).Length() > maxradius)
                maxradius = ( P - A ).Length();

            if ( ( P - B ).Length() > maxradius )
                maxradius = ( P - B ).Length();

            if ( ( P - C).Length() > maxradius )
                maxradius = ( P - C ).Length();
            return maxradius;

        }
        
        /// <summary>
		/// Get the radius of a circumscribed circle to a triangle 
        /// </summary>
        /// <returns></returns>
        public float GetCircleRadius( Vector3 A, Vector3 B, Vector3 C )
        {
            Vector3 AB = ( B - A );
            Vector3 BC = ( C - B );
            Vector3 CA = ( A - C );

            float bclength = BC.Length();
            BC.Normalize();
            AB.Normalize();

            float dotResult = Vector3.Dot( BC, -AB );

            float angleB    = ( float )Math.Acos( Vector3.Dot( BC, -AB ) );

            float value     = bclength / ( 2.0f * ( float )Math.Sin( angleB ) );

            return value;
        }

        public bool IsPointInCircle( Vector3 point, Vector3 circle, float radius )
        {
            return ( ( circle - point ).Length() < (radius-0.01f) );
        }

        public List<Entity> CircleCenter = new List<Entity>();

        /// <summary>
		/// Check an edge and determine if it needs to  be change? If yes the function call itself again
        /// </summary>
        /// <param name="e"></param>
        private void InspectEdge(List<EdgeWE> edges, WingedEdgeMesh mesh)
        {
            foreach(EdgeWE e in edges)
            {
				// Wefirstcheck that the edge connect 2 triangles
                if ( e.LeftFace != null && e.RightFace != null )
                {
                    List<VertexWE> v1 = mesh.GetFaceVertices(e.LeftFace);
                    List<VertexWE> v2 = mesh.GetFaceVertices(e.RightFace);

					//  We get the opposite side of  the left face
                    VertexWE leftOppositeVertex = null;

                    for(int i=0; i< v1.Count; i++)
                    {
                        if( (v1[i] != e.Vertex1) && (v1[i]!=e.Vertex2))
                        {
                            leftOppositeVertex = v1[i];
                        }
                    }

					//We get the opposite side of the right face 
                    VertexWE rightOppositeVertex = null;

                    for(int i=0; i< v2.Count; i++)
                    {
                        if( (v2[i] != e.Vertex1) && (v2[i]!=e.Vertex2))
                        {
                            rightOppositeVertex = v1[i];
                        }
                    }

					// We get back informations about the circles from the left and right face 

                    Vector3 leftCircleCenter = GetCircleCenter( v1[0].Position, v1[1].Position, v1[2].Position );
                    float leftCircleRadius   = GetCircleRadius( v1[0].Position, v1[1].Position, v1[2].Position, leftCircleCenter );

                    Vector3 rightCircleCenter   = GetCircleCenter( v2[0].Position, v2[1].Position, v2[2].Position );
                    float rightCircleRadius     = GetCircleRadius( v2[0].Position, v2[1].Position, v2[2].Position, rightCircleCenter );

					// We check if one of the opposite side can fit in  the circle of  the opposite face 
                    if ( IsPointInCircle( leftOppositeVertex.Position, rightCircleCenter, rightCircleRadius ) || IsPointInCircle( rightOppositeVertex.Position, leftCircleCenter, leftCircleRadius) )
                    {

                        // We switchover the edge
                        List<FaceWE> newfaces = mesh.FlipEdge( e );

                        // We inspect the new edges

                        List<EdgeWE> newEdges = new List<EdgeWE>();

                        for(int i=0; i< newfaces.Count; i++)
                        {
                            for(int j=0; j< newfaces[i].Edges.Count; j++)
                            {
                                if(newEdges.Contains(newfaces[i].Edges[j]))
                                {
                                    newEdges.Add(newfaces[i].Edges[j]);
                                }
                            } 
                        }
                        InspectEdge(newEdges,mesh);
                        return;
                    }
                }
            }
            return ;
        }

        public void DelaunayIncremental()
        {
            List<Entity> p = new List<Entity>();
            List<Entity> currentPoints = new List<Entity>();

            p.AddRange( Points );

            WingedEdgeMesh mesh = new WingedEdgeMesh();

            // We start by making a first triangle

            mesh.AddVertex(-2.0f,-1.0f, 0.0f);
            mesh.AddVertex(-2.0f, 5.0f, 0.0f);
            mesh.AddVertex(2.0f, -1.0f, 0.0f);

            mesh.AddFace( mesh.Vertices[0], mesh.Vertices[1], mesh.Vertices[2] );

			//  Now, we're  going to triangulate vertex per vertex
			// We get back 1 vertex from the set of vertices, we look for the triangle containing it, then,
			// we subdivide  it in 3 triangles. We then check that every new triangle  is Delaunay. If not, we  flip/switchover the
			// opposite edge of  the new vertex
            while ( p.Count > 0 )
            {
                // We start by extracting a point from the set of vertices

                Entity point = p[0];
                p.RemoveAt(0);

                FaceWE f=null;
                int findex = 0;

				// Looking for the triangle that contains the dot
                for(int i=0; i< mesh.Faces.Count; i++)
                {
                    FaceWE currentFace = mesh.Faces[i];
                    List<VertexWE> vertices = mesh.GetFaceVertices(currentFace);

                    bool isIn = Troll3D.TRaycast.PointInTriangle( vertices[0].Position, vertices[1].Position, vertices[2].Position, point.transform_.position_ );

                    if ( isIn )
                    {
                        f=currentFace;
                        findex = i;
                    }
                }

				//  In theory, we should have get back a triangle that contains the vertex. We're now going to divide that 
				// in 3 triangles using the current vertex

                List<FaceWE> faces = mesh.AddVertex( f, point.transform_.position_.X, point.transform_.position_.Y, point.transform_.position_.Z );
                
				// Now we must check  that the triangles we found are Delaunay. Meaning that the 
				// circumscribed circle of the triangle must contains only points of the triangle.
                
                // Checking the faces created
                for ( int i = 0; i < faces.Count; i++ )
                {   
                    InspectEdge(faces[i].Edges,mesh);
                }
            }

			//Removing mesh thing used  for triangulation
            List<FaceWE> faceToRemove = new List<FaceWE>();

            for ( int i = mesh.Faces.Count - 1; i >= 0; i-- )
            {
                if ( mesh.IsFaceBorder( mesh.Faces[i] ) )
                {
                    faceToRemove.Add( mesh.Faces[i] );
                }
            }

            for ( int i = 0; i < faceToRemove.Count; i++ )
            {
                mesh.RemoveFace( faceToRemove[i]);
            }

                for ( int i = 0; i < CircleCenter.Count; i++ )
                {
                    Scene.CurrentScene.RemoveRenderable( ( MeshRenderer )CircleCenter[i].GetComponent( ComponentType.MeshRenderer ) );
                }

                if ( DisplayCenters )
                {
                    for ( int i = 0; i < Circles.Count; i++ )
                    {
                        Scene.CurrentScene.RemoveRenderable( ( LineRenderer )Circles[i].GetComponent( ComponentType.LineRenderer ) );
                    }

                    CircleCenter.Clear();

                    for ( int i = 0; i < mesh.Faces.Count; i++ )
                    {
                        List<VertexWE> vertices = mesh.GetFaceVertices( mesh.Faces[i] );

                        Vector3 circleCenter = GetCircleCenter( vertices[0].Position, vertices[1].Position, vertices[2].Position );
                        float radius = GetCircleRadius( vertices[0].Position, vertices[1].Position, vertices[2].Position, circleCenter );

                        AddCircleCenter( circleCenter.X, circleCenter.Y, radius );
                    }
                }
                else
                {
                    for ( int i = 0; i < Circles.Count; i++ )
                    {
                        Scene.CurrentScene.RemoveRenderable( ( LineRenderer )Circles[i].GetComponent( ComponentType.LineRenderer ) );
                    }
                    CircleCenter.Clear();
                }

            if ( mesh.MakeMesh() != null )
            {
                meshRenderer.model_ = mesh.MakeMesh();
                meshRenderer.SetFillMode( SharpDX.Direct3D11.FillMode.Wireframe );
                meshRenderer.material_.SetMainColor( 1.0f, 0.0f, 0.0f, 1.0f );
            }
            else
            {
            }

            if ( DisplayVoronoi )
            {
                BuildVoronoi( mesh );
            }
            else
            {
                CleanVoronoi();
            }
        }


        public Entity pointDragged = null;
        public bool IsDragging = false;

        public override void OnMouseMove( MouseEvent e )
        {
            if ( IsDragging )
            {
                pointDragged.transform_.position_ = new Vector3( e.mouse_.x / ( float )Screen.Instance.Width * 2 - 1.0f,
                    1.0f - e.mouse_.y / ( float )Screen.Instance.Height * 2, 0.0f );
                DelaunayIncremental();
            }
        }

        public override void OnMouseUp( MouseEvent e )
        {
           
            if ( IsDragging )
            {
                IsDragging = false;
            }
        }

        public override void OnMouseDown( MouseEvent e )
        {
            if ( e.mouse_.leftbutton )
            {
                bool founded = false;
                for ( int i = 0; i < Points.Count && founded==false; i++ )
                {
                    if((Points[i].transform_.position_ - new Vector3(e.mouse_.x / ( float )Screen.Instance.Width * 2 - 1.0f,
                    1.0f - e.mouse_.y / ( float )Screen.Instance.Height * 2, 0.0f )).Length() < 0.01f )
                    {
                        IsDragging = true;
                        pointDragged = Points[i];
                        founded = true;

                    }
                }
            }

            if ( e.mouse_.rightbutton )
            {
                AddPoint( e.mouse_.x / ( float )Screen.Instance.Width * 2 - 1.0f,
                    1.0f - e.mouse_.y / ( float )Screen.Instance.Height * 2 );

                DelaunayIncremental();
            }
        }

        public override void OnKeyDown( KeyboardEvent e )
        {
            if ( e.keycode_ == KeyCode.Key_1 )
            {
                DisplayCenters = !DisplayCenters;
                DelaunayIncremental();
            }
            if ( e.keycode_ == KeyCode.Key_2 )
            {
                DisplayVoronoi = !DisplayVoronoi;
                DelaunayIncremental();
            }
        }
        public void AddCircleCenter( float x, float y, float radius )
        {
            Troll3D.Entity entityCircle = new Troll3D.Entity( Entity );
            entityCircle.transform_.RotateEuler( 0.0f, 3.1415f / 2.0f, 0.0f );
            entityCircle.transform_.Translate( x, y, 0.0f );
            entityCircle.transform_.SetScale( radius, 1.0f, radius );
            LineRenderer lineRenderer = entityCircle.AddComponent<LineRenderer>();


            lineRenderer.material_ = new MaterialDX11( "vDefault.cso", "pUnlit.cso", "gDefaultLine.cso" );
            lineRenderer.material_.SetMainColor( 0.0f, 1.0f, 0.0f, 1.0f );

            lineRenderer.Vertices = Circle.GetLines( 30 );
            lineRenderer.UpdateRenderer();

            Circles.Add( entityCircle );

            Troll3D.Entity entity = new Troll3D.Entity( Entity );

            entity.transform_.SetPosition(
                x,
                y, 
                0.0f );

            entity.transform_.SetScale( 0.02f, 0.02f, 1.0f );

            MaterialDX11 material = new MaterialDX11();
            material.SetMainColor( 0.0f, 1.0f, 0.0f, 1.0f );
            MeshRenderer meshrenderer = entity.AddComponent<MeshRenderer>();
            meshrenderer.material_ = material;
            meshrenderer.model_ = Quad.GetMesh();

            CircleCenter.Add( entity );
        }

        public void AddPoint( float x, float y)
        {
            Troll3D.Entity entity = new Troll3D.Entity( Entity );

            entity.transform_.SetPosition(
                x,
                y,
                0.0f );

            entity.transform_.SetScale( 0.02f, 0.02f, 1.0f );

            MaterialDX11 material = new MaterialDX11();
            material.SetMainColor( 1.0f, 0.0f, 1.0f, 1.0f );
            MeshRenderer meshrenderer = entity.AddComponent<MeshRenderer>();
            meshrenderer.material_ = material;
            meshrenderer.model_ = Quad.GetMesh();

            Points.Add( entity );
        }

        /// <summary>
        /// Contient les points 
        /// </summary>
        public List<Entity> Points = new List<Entity>();
        public List<Entity> Circles = new List<Entity>();

        /// <summary>
        /// Nombre de points à triangulariser
        /// </summary>
        public int PointCount{get;private set;}
        public bool m_start = false;

        public bool DisplayVoronoi=false;
        public bool DisplayCenters = false;
    }
}
