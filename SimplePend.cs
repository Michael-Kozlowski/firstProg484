//============================================================================
// SimplePend.cs  Defines a class for simulating a simple pendulum 
//============================================================================
using System;

namespace Sim 
{
    public class SimplePend
    {
        private double len = 1.1; // pendulum length
        private double g = 9.81;  // gravitation field stregth 
        int n = 2;                // number of states
        private double[] x;       // array of states
        private double[] f;       // right side of the equation value
        //--------------------------------------------------------------------
        // constructor
        //--------------------------------------------------------------------
        public SimplePend()
        {
            x = new double[n];  // array for initial conditions 
            f = new double[n];  // dtheta/dt and dw/dt 

            x[0] = 1.0;  // initial condition for theta in rad
            x[1] = 0.0;  // initial condition for w (omega) in rad/s 
           
        }

        //--------------------------------------------------------------------
        //  step: perform one integration step via Euler's Method
        //--------------------------------------------------------------------
        public void step(double dt)
        {
            rhsFunc(x, f);
            int i; 
            for(i=0;i<n;++i)
            {
                x[i] = x[i] + f[i] * dt;
            }
            //Console.WriteLine($"{f[0].ToString()}   {f[1].ToString()}");
        }
        
        //--------------------------------------------------------------------
        // rhsFunc: Calculate the right handside of the pendulum equations 
        //--------------------------------------------------------------------
        public void rhsFunc(double[] st, double[] ff)
        {
            ff[0] = st[1];
            ff[1] = -g/len * Math.Sin(st[0]);
        }
        //--------------------------------------------------------------------
        // rkinternal: defines the variables w[] and theta[]
        //--------------------------------------------------------------------


        //--------------------------------------------------------------------
        // Getters and Setters
        //--------------------------------------------------------------------
        public double L
        {
            get {return(len);}

            set 
            {
                if(value > 0.0)
                    len = value;
            }
        }

        public double G
        {
            get {return(g);}

            set
            {
                if (value >= 0.0)
                    g = value;
            }
        }

        public double theta
        {
            get{return x[0];}

            set {x[0] = value;}
        }

        public double thetaDot
        {
            get{return x[1];}

            set {x[1] = value;}
        }
    }
}