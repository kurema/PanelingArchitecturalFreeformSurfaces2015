using Rhino;
using Rhino.Geometry;
using Rhino.DocObjects;
using Rhino.Collections;

using System;
using System.IO;
using System.Xml;
using System.Xml.Linq;
using System.Linq;
using System.Data;
using System.Drawing;
using System.Reflection;
using System.Collections;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace RhinoCommonBasic
{
    public class Panelling
    {
        public static object Test() {
            //return new Molds.CubicPolynomial(10, 10, 0.5, 0.4, 0.2, 0.5, 0.2, 0.1).GetSurface();
            var target= Functions.GetSurfaceFromPoints((u, v) => new Point3d(u * 10, v * 10, Math.Sin(u) + Math.Cos(v)), 15, 15);
            var status = new PanellingStatus(target);
            Mold mold1, mold2;
            status.AddMold(mold1=new Molds.Plane(15, 15));
            status.AddMold(mold2=new Molds.CubicPolynomial(15, 15, 0, 0, 0, 0, 0, 0));
            status.AddPanel(new RigidTransformation(-7.5, -7.5, 0, 0, 0, 0), mold1);
            status.AddPanel(new RigidTransformation(-7.5, 7.5, 0, 0, 0, 0), mold2);
            status.AddPanel(new RigidTransformation(7.5, -7.5, 0, 0, 0, 0), mold2);
            status.AddPanel(new RigidTransformation(7.5, 7.5, 0, 0, 0, 0), mold2);
            status.CurveNetworks.AddRandomPointsOnSurface(target, 10);
            status.TryOptimizationTransform(1);
            return status.GetPanels();
        }

        public interface Mold : ParamaterProvider
        {
            Surface GetSurface();
        }

        public interface ParamaterProvider
        {
            IParamaters Paramaters { get; }
        }

        public interface TestProvider
        {
            double[] GetTest();
        }

        public interface JacobianProvider
        {
            Matrix GetJacobian();
        }

        public class PanellingStatus: TestProvider, JacobianProvider,ParamaterProvider
        {
            public Surface TargetSurface;
            public CurveNetworks CurveNetworks = new CurveNetworks();
            public Mold[] MoldDepot { get { return _MoldDepot; } set { _MoldDepot = value; InitParamaters(); } }
            private Mold[] _MoldDepot=new Mold[0];
            public RigidTransformation[] Transformations { get; private set; }
            public Mold[] PanelMoldSelection { get; private set; }

            public PanellingStatus(Surface Target)
            {
                TargetSurface = Target;
                Transformations = new RigidTransformation[0];
                PanelMoldSelection = new Mold[0];

                var para = new ParamatersCombination();
                para.Init(CurveNetworks.Paramaters);
                Paramaters = para;
            }

            public void InitParamaters()
            {
                var result= new List<IParamaters>();
                result.Add(CurveNetworks.Paramaters);
                foreach (var item in MoldDepot)
                {
                    result.Add(item.Paramaters);
                }
                foreach(var item in Transformations)
                {
                    result.Add(item.Paramaters);
                }

                var para = new ParamatersCombination();
                para.Init(result.ToArray());
                Paramaters = para;

            }

            public Surface GetPanel(int i)
            {
                var surface = PanelMoldSelection[i].GetSurface();
                surface.Transform(Transformations[i].ToTransform());
                return surface;
            }

            public Surface[] GetPanels()
            {
                Surface[] result = new Surface[Transformations.Count()];
                for(int i = 0; i < Transformations.Count(); i++)
                {
                    result[i]=GetPanel(i);
                }
                return result;
            }

            public void AddMold(Mold item)
            {
                var temp = MoldDepot.ToList();
                temp.Add(item);
                MoldDepot = temp.ToArray();

                InitParamaters();
            }

            public void RemoveMold(Mold item)
            {
                var temp = MoldDepot.ToList();
                var idx= temp.IndexOf(item);
                temp.RemoveAt(idx);
                MoldDepot = temp.ToArray();

                InitParamaters();
            }

            public void AddPanel(RigidTransformation tr,Mold mold)
            {
                if (!MoldDepot.Contains(mold))
                {
                    AddMold(mold);
                }
                var ttr = Transformations.ToList();
                ttr.Add(tr);
                Transformations = ttr.ToArray();

                var tpm = PanelMoldSelection.ToList();
                tpm.Add(mold);
                PanelMoldSelection = tpm.ToArray();

                InitParamaters();
            }

            public IParamaters Paramaters
            {
                get;private set;
            }

            public double[] GetTest()
            {
                var result = new List<double>();
                //result.AddRange(TestSurfaceFitting());
                result.AddRange(TestDivergence());
                return result.ToArray();
            }

            public Matrix GetJacobian()
            {
                return GetJacobian(Paramaters);
            }

            public Matrix GetJacobian(IParamaters p)
            {
                var orgtest = GetTest();
                Matrix result = new Matrix(orgtest.Count(),p.Count());
                for(int i = 0; i < p.Count(); i++)
                {
                    double dx = 1e-6;
                    double val= p.Get(i);
                    val += dx;//微分時のdx
                    p.Set(i, val);
                    var newtest= GetTest();
                    val -= dx;
                    for(int j = 0; j < orgtest.Count(); j++)
                    {
                        result[j, i] = (newtest[j]-orgtest[j]) / dx;
                    }
                }
                return result;
            }

            public double[] TryOptimization(int count = 1,double Scale=1.0)
            {
                List<double> history = new List<double>();

                Matrix BestMatrix = Functions.ArrayToMatrix(Paramaters.ToArray());//不要
                //double BestValue = double.MaxValue;
                for (int i = 0; i < count; i++)
                {
                    var paraMatrix = Functions.ArrayToMatrix(Paramaters.ToArray());
                    var test = GetTest();
                    Matrix m = Functions.GaussNewtonMethod(GetJacobian(), Functions.ArrayToMatrix(test), paraMatrix);
                    //m.Scale(Scale);
                    //Paramaters.Init(Functions.MatrixToArray(m + paraMatrix));
                    TryLengthEstimation( Functions.MatrixToArray(paraMatrix), Functions.MatrixToArray(m));

                    var sqsum = Functions.GetSquareSum(test);
                    //if (sqsum < BestValue) { BestMatrix = m + paraMatrix; }
                    history.Add(sqsum);
                }
                //Paramaters.Init(Functions.MatrixToArray(BestMatrix));

                return history.ToArray();
            }

            public void TryLengthEstimation(double[] org,double[] deg)
            {
                IParamaters[] paras = Paramaters.GetChildParamaters();
                int cnt = 0;

                int searchWidthL = 5;
                int searchWidthN = 3;
                int searchCountN = 5;
                foreach (var item in paras)
                {
                    double bestScale = 0;
                    double bestValue = double.MaxValue;
                    double rangeMax = 1e5;
                    double rangeMin = 1e-5;
                    double rangeMaxNew = 1e5;
                    double rangeMinNew = 1e-5;

                    double searchStep = Math.Pow(10, Math.Log10(rangeMax / rangeMin) / searchWidthL);

                        for (double val = rangeMin; val < rangeMax; val *= searchStep)
                        {
                            double scale = val;

                            TryLengthEstimationTestScale(item, scale, cnt, org, deg);
                            var sqsum = Functions.GetSquareSum(GetTest());
                            if (sqsum < bestValue)
                            {
                                bestValue = sqsum;
                                bestScale = scale;
                                rangeMaxNew = scale * Math.Sqrt(searchStep);
                                rangeMinNew = scale / Math.Sqrt(searchStep);
                            }
                        }
                    rangeMax = rangeMaxNew;
                    rangeMin = rangeMinNew;

                    for (int i = 0; i < searchCountN; i++) {
                        for (int j = 0; j <= searchWidthN; j++)
                        {
                            double scale = rangeMin + (rangeMax - rangeMin) / searchWidthN * j;

                            TryLengthEstimationTestScale(item, scale, cnt, org, deg);
                            var sqsum = Functions.GetSquareSum(GetTest());
                            if (sqsum < bestValue)
                            {
                                bestValue = sqsum;
                                bestScale = scale;
                                rangeMaxNew = scale + (rangeMax - rangeMin) / searchWidthN / 2.0;
                                rangeMinNew = scale - (rangeMax - rangeMin) / searchWidthN / 2.0;
                            }
                        }
                        
                        rangeMax = rangeMaxNew;
                        rangeMin = rangeMinNew;
                    }
                    TryLengthEstimationTestScale(item, bestScale, cnt, org, deg);

                    cnt += item.Count();
                }
            }

            protected void TryLengthEstimationTestScale(IParamaters item,double Scale,int cnt,double[] org,double[] deg)
            {
                for (int i = 0; i < item.Count(); i++)
                {
                    item.Set(i, org[i+cnt]+Scale*deg[i+cnt]);
                }
            }

            public double[] TryOptimizationTransform(int count = 1, double Scale = 1.0)
            {
                List<double> history = new List<double>();

                for (int i = 0; i < count; i++)
                {
                    var test = GetTest();
                    for (int j = 0; j < Transformations.Count(); j++)
                    {
                        var paraMatrix = Functions.ArrayToMatrix(this.Transformations[j].Paramaters.ToArray());
                        Matrix m = Functions.GaussNewtonMethod(GetJacobian(Transformations[j].Paramaters), Functions.ArrayToMatrix(test), paraMatrix);
                        m.Scale(Scale);
                        this.Transformations[j].Paramaters.Init(Functions.MatrixToArray(m + paraMatrix));
                    }
                    history.Add(Functions.GetSquareSum(test));
                }
                return history.ToArray();
            }


            public double[] TestSurfaceFitting()
            {
                List<double> result = new List<double>();
                var controlPoints = CurveNetworks.ControlPoints;
                foreach(var cp in controlPoints)
                {
                    result.Add(Functions.GetPointToSurfaceDistance(cp,TargetSurface));
                }
                return result.ToArray();
            }

            public double[] TestCost()
            {
                List<double> result = new List<double>();
                foreach(var item in this.Transformations)
                {
                    result.AddRange(item.Paramaters.ToArray());
                }
                return result.ToArray();
            }

            public double[] TestDivergence(params Surface[] panels) {
                //fix me
                var panellist = panels.Count() == 0 ? GetPanels() : panels;
                var controlPoints = CurveNetworks.ControlPoints;
                double Length1st = double.MaxValue;
                double Length2nd = double.MaxValue;
                List<double> result = new List<double>();

                foreach (var cp in controlPoints)
                {
                    foreach(var panel in panellist)
                    {
                        double dist = Functions.GetPointToSurfaceDistance(cp, panel);
                        if (Length1st > dist)
                        {
                            Length2nd = Length1st;
                            Length1st = dist;
                        }
                        else if(Length2nd > dist){
                            Length2nd = dist;
                        }
                    }
                    result.Add(Length1st);
                    result.Add(Length2nd);
                }
                return result.ToArray();
            }
        }

        public class CurveNetworks: ParamaterProvider
        {
            public Point3d[] ControlPoints { get
                {
                    List<Point3d> result = new List<Point3d>();
                    for(int i = 0; i < Paramaters.Count(); i += 3)
                    {
                        result.Add(new Point3d(Paramaters.Get(i), Paramaters.Get(i+1), Paramaters.Get(i+2)));
                    }
                    return result.ToArray();
                }
            }

            public IParamaters Paramaters
            {
                get; private set;
            }

            public CurveNetworks()
            {
                Paramaters = new Paramaters(0);
            }

            public void AddPoint(params Point3d[] ps)
            {
                var org= Paramaters.ToArray().ToList();
                foreach (var p in ps)
                {
                    org.Add(p.X);
                    org.Add(p.Y);
                    org.Add(p.Z);
                }
                Paramaters.Init(org.ToArray());
            }

            public void AddRandomPointsOnSurface(Surface surface,int count)
            {
                AddPoint(Functions.GetRandomPointsOnSurface(surface, count));
            }

            public void InitControlPoints(Curve[] arg, int Count)
            {
                var result = new List<double>();
                foreach (Curve item in arg)
                {
                    for (int i = 0; i < Count; i++) {
                        var point= (item.PointAtNormalizedLength(i/(Count-1)));
                        result.Add(point.X);
                        result.Add(point.Y);
                        result.Add(point.Z);
                    }
                }
                Paramaters = new Paramaters(result.Count()); ;
                Paramaters.Init(result.ToArray());
            }
        }

        public interface IParamaters : IEnumerable
        {
            event EventHandler ValueChanged;
            void Init(params double[] arg);
            bool Set(int target, double value);
            double[] ToArray();
            double Get(int target);
            int Count();
            IParamaters[] GetChildParamaters();
        }

        public class ParamatersCombination : IParamaters
        {
            public event EventHandler ValueChanged;
            public IParamaters[] ParamatersMember { get; private set; }

            protected virtual void OnValueChanged(EventArgs e)
            {
                if (ValueChanged != null) ValueChanged(this, e);
            }

            public double Get(int target)
            {
                int cnt = 0;
                foreach (var item in ParamatersMember)
                {
                    if (cnt <= target && target < cnt + item.Count())
                    {
                        return item.Get(target - cnt);
                    }
                    cnt += item.Count();
                }
                throw new Exception();
            }

            public IEnumerator GetEnumerator()
            {
                foreach (var item in ParamatersMember)
                {
                    var itemAr= item.ToArray();
                    foreach(var item2 in itemAr)
                    {
                        yield return item2;
                    }
                }
            }

            public void Init(params double[] arg)
            {
                int cnt = 0;
                foreach (var item in this.GetChildParamaters())
                {
                    for(int i = 0; i < item.Count(); i++)
                    {
                        item.Set(i,arg[cnt]);
                        cnt++;
                    }
                }
                //長すぎる配列分は無視。
            }

            public void Init(params IParamaters[] Ps)
            {
                ParamatersMember = Ps;
            }

            public void AddMember(params Paramaters[] Ps)
            {
                var para= ParamatersMember.ToList();
                para.AddRange(Ps);
                ParamatersMember = para.ToArray();
            }

            public bool Set(int target, double value)
            {
                int cnt = 0;
                foreach (var item in ParamatersMember)
                {
                    if (cnt <= target && target < cnt + item.Count())
                    {
                        item.Set(target - cnt, value);
                        return true;
                    }
                    cnt += item.Count();
                }
                return false;
            }

            public double[] ToArray()
            {
                List<double> Result = new List<double>();
                foreach (var item in ParamatersMember)
                {
                    Result.AddRange(item.ToArray());
                }
                return Result.ToArray();
            }

            public int Count()
            {
                int result = 0;
                foreach(var item in ParamatersMember)
                {
                    result += item.Count();
                }
                return result;
            }

            public IParamaters[] GetChildParamaters()
            {
                List<IParamaters> result = new List<IParamaters>();
                foreach(var item in this.ParamatersMember)
                {
                    result.AddRange(item.GetChildParamaters());
                }
                return result.ToArray();
            }
        }

        public class Paramaters : IParamaters
        {
            protected double[] Content;

            public event EventHandler ValueChanged;

            protected virtual void OnValueChanged(EventArgs e)
            {
                if (ValueChanged != null) ValueChanged(this, e);
            }

            public Paramaters(int Count)
            {
                Content = new double[Count];
            }

            public IEnumerator GetEnumerator()
            {
                foreach(var item in Content)
                {
                    yield return item;
                }
            }

            public void Init(params double[] arg)
            {
                Content = arg;
            }

            public bool Set(int target,double value)
            {
                if (Count() > target)
                {
                    Content[target] = value;
                    return true;
                }
                else
                {
                    return false;
                }
            }

            public double[] ToArray()
            {
                return Content;
            }

            public double Get(int target)
            {
                return Content[target];
            }
            public int Count()
            {
                return Content.Count();
            }

            public IParamaters[] GetChildParamaters()
            {
                return new IParamaters[] { this };
            }
        }

        public class Molds
        {
            public abstract class PointSurfaceBase : Mold
            {
                public abstract IParamaters Paramaters { get; }

                public int uCnt = 10, vCnt = 10;
                public Surface GetSurface() {
                    return Functions.GetSurfaceFromPoints((u, v) => { return GetPoint(u, v); },uCnt,vCnt);
                }
                public abstract Point3d GetPoint(double u, double v);
            }

            public class Plane : Mold
            {
                public IParamaters Paramaters
                {
                    get;private set;
                }

                public Plane(double Width,double Height)
                {
                    Paramaters = new Paramaters(2);
                    Paramaters.Set(0, Width);
                    Paramaters.Set(1, Height);
                }

                public Surface GetSurface()
                {
                    return new PlaneSurface(Rhino.Geometry.Plane.WorldXY, new Interval(-Paramaters.Get(0) / 2.0, Paramaters.Get(0) / 2.0), new Interval(-Paramaters.Get(1) / 2.0, Paramaters.Get(1) / 2.0));
                }
            }

            public class CubicPolynomial : PointSurfaceBase
            {
                public CubicPolynomial(double width,double height,double a,double b,double c,double d,double e,double f)
                {
                    Paramaters.Init(width, height, a, b, c, d, e, f);
                }

                public override IParamaters Paramaters
                {
                    get { return _param; }
                }
                private Paramaters _param = new Paramaters(8);
 
                public override Point3d GetPoint(double u, double v)
                {
                    double[] p = Paramaters.ToArray();
                    return new Point3d(u * p[0] / 2.0, v * p[1] / 2.0, (p[2] * u * u + p[3] * v * v + p[4] * u * u * u + p[5] * u * u * v + p[6] * u * v * v + p[7] * v * v * v) * (p[0] + p[1]) / 2.0);
                }

                public double[] Get6DParamater(int Count)
                {
                    double[] p = Paramaters.ToArray();
                    return new double[]
                    {
                        2/3/Math.Sqrt(5)*p[2],2/3/Math.Sqrt(5)*p[3],(p[5]+p[7])/Math.Sqrt(15),(p[4]+p[6])/Math.Sqrt(15),Math.Sqrt(8.0/15.0)*p[4],Math.Sqrt(8.0/15.0)*p[7]
                    };
                }
            }
        }

        public class RigidTransformation : ParamaterProvider
        {
            //x,y,z,rx,ry,rz
            public IParamaters Paramaters
            {
                get;private set;
            }

            public RigidTransformation(double x,double y,double z,double rx,double ry,double rz)
            {
                var pc = new ParamatersCombination();
                pc.Init(new Paramaters(3), new Paramaters(3));
                pc.Init(x, y, z, rx, ry, rz);
                Paramaters = pc;
            }

            public Transform ToTransform()
            {
                double[] p = Paramaters.ToArray();
                return Transform.Multiply( Transform.Translation(p[0], p[1], p[2])
                    , Transform.Multiply(Transform.Rotation(Math.Atan(p[5]), Vector3d.ZAxis, Point3d.Origin)
                    , Transform.Multiply(Transform.Rotation(Math.Atan(p[4]), Vector3d.YAxis, Point3d.Origin)
                    , Transform.Rotation(Math.Atan(p[3]), Vector3d.XAxis, Point3d.Origin))));
            }

        }

        public static class Functions
        {
            public static Surface GetSurfaceFromPoints(Func<double,double,Point3d> func,int uCnt,int vCnt)
            {
                List<Point3d> points = new List<Point3d>();
                for (int i = 0; i < uCnt;i++) {
                    for (int j = 0; j < vCnt; j++)
                    {
                        points.Add(func((double)i / (uCnt - 1) * 2 - 1, (double)j / (vCnt - 1) * 2 - 1));
                    }
                }
                return (Surface)NurbsSurface.CreateFromPoints(points.ToArray(), uCnt, vCnt, 3, 3);
            }

            static public double GetPointToSurfaceDistance(Point3d point,Surface surface)
            {
                double u, v;
                surface.ClosestPoint(point, out u, out v);
                return surface.PointAt(u, v).DistanceTo(point);
            }

            static public Matrix ArrayToMatrix(double[] arg)
            {
                var result = new Matrix(arg.Count(),1);
                for(int i = 0; i < arg.Count(); i++)
                {
                    result[i, 0] = arg[i];
                }
                return result;
            }

            static public double[] MatrixToArray(Matrix arg)
            {
                List<double> result = new List<double>();
                for (int j = 0; j < arg.RowCount; j++)
                {
                    for (int i = 0; i < arg.ColumnCount; i++)
                    {
                        result.Add(arg[j, i]);
                    }
                }
                return result.ToArray();
            }

            static public Point3d[] GetRandomPointsOnSurface(Surface surface, int count) {
                Random rd = new Random();
                var result = new Point3d[count];
                for(int i=0;i< count; i++)
                {
                    var dom1 = surface.Domain(0);
                    var dom2 = surface.Domain(1);
                    result[i] = surface.PointAt(dom1.Min + dom1.Length * rd.NextDouble(), dom2.Min + dom2.Length * rd.NextDouble());
                    
                }
                return result;
            }

            static public double GetSquareSum(double[] arg)
            {
                double result = 0;
                foreach(var item in arg) { result += item;  }
                return result;
            }

            static public Polyline GetSimpleGraph(double[] d)
            {
                double maxValue=double.MinValue;
                double minValue=double.MaxValue;
                foreach(var item in d)
                {
                    maxValue = Math.Max(maxValue, item);
                    minValue = Math.Min(minValue, item);
                }
                double lenValue = maxValue - minValue;
                var pl = new Polyline(d.Count());
                for (int i = 0; i < d.Count(); i++)
                {
                    //pl.Add(new Point3d((double)i / (d.Count() - 1.0), (d[i] - minValue) / lenValue, 0));
                    pl.Add(new Point3d((double)i / (d.Count() - 1.0), d[i]/maxValue, 0));
                }

                return pl;
            }

            static public Matrix GaussNewtonMethod(Matrix Jacobian, Matrix EnergyVector, Matrix StatusVector, double ZeroTolerance = 1e-10)
            {
                Matrix M1 = Jacobian.Duplicate();
                M1.Transpose();
                Matrix M2 = M1 * Jacobian;
                M2.Invert(ZeroTolerance);
                M1 = M1 * EnergyVector;
                M1 = M2 * M1;
                M1.Scale(-1);
                return M1;
            }
        }
    }
}
