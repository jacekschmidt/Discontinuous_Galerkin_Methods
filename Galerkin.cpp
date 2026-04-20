#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>

// advection element class
class AdvectionElement
{
    public:
        // Pointer to the left neighbour
        AdvectionElement *Left_neighbour_pt;

        // Pointer to the right neighbour
        AdvectionElement *Right_neighbour_pt;
        
        // Storage for the coordinates
        MVector X;

        // Storage for the unknowns
        MVector U;

        MVector Utemp;

        MVector K1;

        MVector K2;

        MVector K3;

        MVector K4;
        
        // Constructor: initialise the vectors to hold two entries.
        AdvectionElement(): Left_neighbour_pt(nullptr), Right_neighbour_pt(nullptr), X(2), U(2), Utemp(2), K1(2), K2(2), K3(2), K4(2)
        {}

        virtual ~AdvectionElement() = default;

        // Return the value of the coordinate at local coordinate s using
        // equation (1.2)
        double interpolated_x(double s) 
        {return 0.5*(X[0]+X[1]+s*(X[1]-X[0]));}

        // Return the value of the unknown at local coordinate s using
        // equation (1.4)
        double interpolated_u(double s) 
        {return 0.5*(Utemp[0]+Utemp[1]+s*(Utemp[1]-Utemp[0]));}

        // Calculate the flux
        virtual double flux(double u)
        {
        return u;
        }

        // flux derivative
        virtual double max_flux_derivative(double a, double b)
        {
        return 1.0;
        }

        // Calculate the integral of the flux function over the element
        // using the two-point Gauss rule
        double integrate_flux()
        {
        return flux(interpolated_u(-1/std::sqrt(3)))+flux(interpolated_u(1/std::sqrt(3)));
        }

        // numerical flux 
        virtual double h(double a, double b)
        {
        return 0.5*(flux(a)+flux(b)-max_flux_derivative(a, b)*(b-a));
        }

        // helpful for rk4
        double F_0_plus_h_0()
        {
            return -0.5*integrate_flux()+h((*Left_neighbour_pt).Utemp[1],Utemp[0]);
        }

        // helpful for rk4
        double F_1_plus_h_1()
        {
            return 0.5*integrate_flux()-h(Utemp[1],(*Right_neighbour_pt).Utemp[0]);
        }

        // udot
        MVector udot(double minbound, double maxbound, int N)
        {
            double stepsize=(maxbound-minbound)/N;
            MVector UD(2);
            UD[0]=(2.0/stepsize)*(2*(F_0_plus_h_0())-(F_1_plus_h_1()));
            UD[1]=(2.0/stepsize)*(-(F_0_plus_h_0())+2*(F_1_plus_h_1()));
            return UD;
        }

        // updating function with Eiger 1
        void timestep(double dt, double minbound, double maxbound, int N)
        {
        U=Utemp+dt*udot(minbound, maxbound, N);
        }

        // updating function with RK4
        void timestep_RK4(double dt)
        {
        U=U+(dt/6.0)*(K1+2*K2+2*K3+K4);
        }

        // minmod function
        double minmod(double sl, double sc, double sr, double a)
        {
            if (sl > 0.0 && sc > 0.0 && sr > 0.0)
            {
                return std::min({a*sl, sc, a*sr});
            }
            if (sl < 0.0 && sc < 0.0 && sr < 0.0)
            {
                return std::max({a*sl, sc, a*sr});
            }
            return 0.0;
        }

        void minmod_limiter(double a)
        {
            double ubar  = 0.5 * (Utemp[0] + Utemp[1]);
            double ubarL = 0.5 * ((*Left_neighbour_pt).Utemp[0] +
                                (*Left_neighbour_pt).Utemp[1]);
            double ubarR = 0.5 * ((*Right_neighbour_pt).Utemp[0] +
                                (*Right_neighbour_pt).Utemp[1]);
            double sl = ubar - ubarL;
            double sr = ubarR - ubar;
            double sc = Utemp[1] - Utemp[0];
            double s_new = minmod(sl, sc, sr, a);
            Utemp[0] = ubar - 0.5 * s_new;
            Utemp[1] = ubar + 0.5 * s_new;
        }


}; //End of the class definition

// Burger's element class
class BurgersElement: public AdvectionElement
{
    public:
        // Calculate the flux
            virtual double flux(double u)
            {
            return 0.5*(u*u);
            }

        // flux derivative
            virtual double max_flux_derivative(double a, double b)
            {
            return std::max(std::abs(a),std::abs(b));
            }
};

// galerkin solver with euler 1
void Galerkin_Solver(std::vector<AdvectionElement*> &elements, double dt, double minbound, double maxbound)
{
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->Utemp=elements[i]->U;
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->timestep(dt, minbound, maxbound, elements.size());
    } 
}

// galerkin solver with euler 1 and minmod
void Galerkin_Solver_minmod(std::vector<AdvectionElement*> &elements, double dt, double minbound, double maxbound, double a)
{
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->Utemp=elements[i]->U;
    }
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->minmod_limiter(a);
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->timestep(dt, minbound, maxbound, elements.size());
    } 
}

// implements 1 Galerkin Iteration with RK4
void Galerkin_Solver_RK4(std::vector<AdvectionElement*> &elements, double dt, double minbound, double maxbound)
{
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->Utemp=elements[i]->U;
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->K1=elements[i]->udot(minbound, maxbound, elements.size());
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->Utemp=elements[i]->U+(dt/2.0)*elements[i]->K1;
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->K2=elements[i]->udot(minbound, maxbound, elements.size());
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->Utemp=elements[i]->U+(dt/2.0)*elements[i]->K2;
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->K3=elements[i]->udot(minbound, maxbound, elements.size());
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->Utemp=elements[i]->U+dt*elements[i]->K3;
    }
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->K4=elements[i]->udot(minbound, maxbound, elements.size());
    } 
    for (int i=0; i<elements.size(); i++)
    {
        elements[i]->timestep_RK4(dt);
    }  
}

// implements 1 Galerkin Iteration with RK4 and minmod limiter
void Galerkin_Solver_RK4_minmod(std::vector<AdvectionElement*> &elements,
                                double dt, double minbound, double maxbound, double a)
{
    // ---- Stage 1 ----
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->Utemp = elements[i]->U;
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->minmod_limiter(a);
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->K1 = elements[i]->udot(minbound, maxbound, elements.size());
    }

    // ---- Stage 2 ----
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->Utemp = elements[i]->U + (dt/2.0) * elements[i]->K1;
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->minmod_limiter(a);
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->K2 = elements[i]->udot(minbound, maxbound, elements.size());
    }

    // ---- Stage 3 ----
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->Utemp = elements[i]->U + (dt/2.0) * elements[i]->K2;
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->minmod_limiter(a);
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->K3 = elements[i]->udot(minbound, maxbound, elements.size());
    }

    // ---- Stage 4 ----
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->Utemp = elements[i]->U + dt * elements[i]->K3;
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->minmod_limiter(a);
    }

    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->K4 = elements[i]->udot(minbound, maxbound, elements.size());
    }

    // ---- Final RK4 update ----
    for (int i = 0; i < elements.size(); i++)
    {
        elements[i]->timestep_RK4(dt);
    }
}

// sin mesh constructor
void sinconstructor(std::vector<AdvectionElement*> &elements, double minbound, double maxbound)
{
    double h=(2.0*M_PI)/elements.size();
    for (int i=0; i<elements.size(); i++)
        {
            elements[i]->X[0] = i*h ;
            elements[i]->X[1] = (i+1)*h;
            elements[i]->U[0] = 1.5 + sin(elements[i]->X[0]);
            elements[i]->U[1] = 1.5 + sin(elements[i]->X[1]);
            if (i > 0)
                elements[i]->Left_neighbour_pt = elements[i-1];
            else
                elements[i]->Left_neighbour_pt = elements[elements.size()-1];  // boundary

            if (i <elements.size()-1)
                elements[i]->Right_neighbour_pt = elements[i+1];
            else
                elements[i]->Right_neighbour_pt = elements[0]; // boundary
        }
}

// square wave constructor
void squareconstructor(std::vector<AdvectionElement*> &elements, double minbound, double maxbound)
{
    double h=2.0/elements.size();
    for (int i=0; i<elements.size(); i++)
        {
            elements[i]->X[0] = i*h ;
            elements[i]->X[1] = (i+1)*h;
            if(i<=elements.size()/2){elements[i]->U[0] = 1.0;}
            else {elements[i]->U[0] = 0.0;}
            if(i+1<=elements.size()/2){elements[i]->U[1] = 1.0;}
            else {elements[i]->U[1] = 0.0;}
            if (i > 0)
                elements[i]->Left_neighbour_pt = elements[i-1];
            else
                elements[i]->Left_neighbour_pt = elements[elements.size()-1];  // boundary

            if (i <elements.size()-1)
                elements[i]->Right_neighbour_pt = elements[i+1];
            else
                elements[i]->Right_neighbour_pt = elements[0]; // boundary
        }
}

// plots the progression of an initialisation
void progression_plotter(const std::string& filename, std::vector<AdvectionElement*> &elements, double T, double dt, double minbound, double maxbound)
{
    std::ofstream demofile;
    demofile.open(filename);
    if (!demofile)
    {
        std::cout<<"unable to open file";
        return;
    }
    int endtime=ceil(T/dt)+1;
    for (int i=0; i<endtime; i++)
    {
        for (int j=0; j<elements.size(); j++)
        {
            demofile << elements[j]->interpolated_u(0.0) << " ";
        }
        Galerkin_Solver(elements, dt, minbound, maxbound);
    }
    demofile.close();
}

// residual computation for sin
double sin_adv_error(std::vector<AdvectionElement*> &elements, double T, double minbound, double maxbound)
{
    double sum=0.0;
    double h=(maxbound-minbound)/elements.size();
    for(int i=0; i<elements.size(); i++)
    {
        sum+=pow(elements[i]->interpolated_u(0.0)-(1.5+sin(elements[i]->interpolated_x(0.0)-T)),2);
    }
    return std::sqrt(h*sum);
} 

int main()
{
    /* testing Advection element class initilisation
    {
        int N = 7;
        double h=(2.0*M_PI)/N;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        for (int i=0; i<N; i++)
        {
            elements[i]->X[0] = i*h ;
            elements[i]->X[1] = (i+1)*h;
            elements[i]->U[0] = 1.5 + sin(elements[i]->X[0]);
            elements[i]->U[1] = 1.5 + sin(elements[i]->X[1]);
            elements[i]->Utemp=elements[i]->U;
            if (i > 0)
                elements[i]->Left_neighbour_pt = elements[i-1];
            else
                elements[i]->Left_neighbour_pt = elements[N-1];  // boundary

            if (i < N-1)
                elements[i]->Right_neighbour_pt = elements[i+1];
            else
                elements[i]->Right_neighbour_pt = elements[0]; // boundary

            // write code to fill out other member data for the i'th element
            // here
        }
        for (int i=0; i<elements.size(); i++)
        {
            std::cout << elements[i]->interpolated_x(0.0) << "," << elements[i]->interpolated_u(0.0) << std::endl;
        }
    }
    */
    
    /* advection with sin
    {
        int N=100;
        double dt=0.001;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2*M_PI;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        sinconstructor(elements,minbound,maxbound);
        for (int i=0; i<1000; i++)
        {
            Galerkin_Solver(elements, dt, minbound, maxbound);
        }
        for (int i=0; i<elements.size(); i++)
        {
            std::cout << elements[i]->interpolated_u(0.0) << std::endl;
        }
    }
    */

    /* progression plotter with sin
    {
        int N=100;
        double dt=0.001;
        double minbound=0.0;
        double maxbound=2*M_PI;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        sinconstructor(elements,minbound,maxbound);
        progression_plotter("progression_of_sin_adv.txt", elements, 1.0, dt, minbound, maxbound);
    }
    */

    /* error dependence on timestep
    {
        int N=500;
        double minbound=0.0;
        double maxbound=2*M_PI;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        sinconstructor(elements,minbound,maxbound);
        for(int i=3750; i>3500; i--)
        {
            double dt=1.0/i;
            for (int j=0; j<i; j++)
            {
                Galerkin_Solver(elements, dt, minbound, maxbound);
            }
            std::cout<<sin_adv_error(elements, 1.0, minbound, maxbound)<<std::endl;
            sinconstructor(elements,minbound,maxbound);
        }        
    }
    */

    /* error dependence on h
    {
        double dt=1.0/10;
        double minbound=0.0;
        double maxbound=2*M_PI;
        for(int i=60; i>59; i--)
        {
            int N=i;
            std::vector<AdvectionElement*> elements(N);
            for (int i = 0; i < N; ++i)
            {
                elements[i] = new AdvectionElement();
            }
            sinconstructor(elements,minbound,maxbound);
            for (int j=0; j<10; j++)
            {
                Galerkin_Solver(elements, dt, minbound, maxbound);
            }
            std::cout<<sin_adv_error(elements, 1.0, minbound, maxbound)<<std::endl;
        }        
    }
    */

    /* error depends on dt and h
    {
        std::ofstream demofile;
        demofile.open("Error_surface.txt");
        if (!demofile)
        {
            std::cout<<"unable to open file";
            return 1;
        }
        double minbound=0.0;
        double maxbound=2*M_PI;
        for(int i=6000; i>=1200; i=i-600)
        {
            int N=i;
            std::vector<AdvectionElement*> elements(N);
            for (int i = 0; i < N; ++i)
            {
                elements[i] = new AdvectionElement();
            }
            sinconstructor(elements,minbound,maxbound);
            for (int j=1000; j>=200; j=j-100)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }
            for (int j=100; j>=30; j=j-10)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }
            for (int j=20; j>=10; j--)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }        
        }
        for(int i=600; i>=180; i=i-60)
        {
            int N=i;
            std::vector<AdvectionElement*> elements(N);
            for (int i = 0; i < N; ++i)
            {
                elements[i] = new AdvectionElement();
            }
            sinconstructor(elements,minbound,maxbound);
            for (int j=1000; j>=200; j=j-100)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }
            for (int j=100; j>=30; j=j-10)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }
            for (int j=20; j>=10; j--)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }        
        }
        for(int i=120; i>=60; i=i-6)
        {
            int N=i;
            std::vector<AdvectionElement*> elements(N);
            for (int i = 0; i < N; ++i)
            {
                elements[i] = new AdvectionElement();
            }
            sinconstructor(elements,minbound,maxbound);
            for (int j=1000; j>=200; j=j-100)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }
            for (int j=100; j>=30; j=j-10)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }
            for (int j=20; j>=10; j--)
            {
                double dt=1.0/j;
                for (int k=0; k<j; k++)
                {
                    Galerkin_Solver(elements, dt, minbound, maxbound);
                }
                double a=sin_adv_error(elements, 1.0, minbound, maxbound);
                if (a>0.4) {demofile<<0.4<<" ";}
                else {demofile<<a<<" ";}
                sinconstructor(elements,minbound,maxbound);
            }            
        }  
        demofile.close();
    }
    */

    /* square constructor test
    {
        int N=3;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        squareconstructor(elements,minbound, maxbound);
        for (int i=0; i<elements.size(); i++)
        {
            std::cout << elements[i]->interpolated_x(0.0) << "," << elements[i]->interpolated_u(0.0) << std::endl;
        }
    }
    */

    /* advection with square
    {
        int N=100;
        double dt=0.001;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        squareconstructor(elements,minbound,maxbound);
        for (int i=0; i<1000; i++)
        {
            Galerkin_Solver(elements, dt, minbound, maxbound);
        }
        for (int i=0; i<elements.size(); i++)
        {
            elements[i]->Utemp=elements[i]->U;
            std::cout << elements[i]->interpolated_u(0.0) <<  std::endl;
        }
    }
    */
    
    
    /* progression with square
    {
        int N=100;
        double dt=0.00001;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        squareconstructor(elements,minbound,maxbound);
        progression_plotter("progression_of_square_adv.txt", elements, 1.0, dt, minbound, maxbound);
    }
    */

    /* advection with square and minmod
    {
        int N=100;
        double dt=0.001;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        squareconstructor(elements,minbound,maxbound);
        for (int i=0; i<1000; i++)
        {
            Galerkin_Solver_RK4_minmod(elements, dt, minbound, maxbound, 2.0);
        }
        for (int i=0; i<elements.size(); i++)
        {
            elements[i]->Utemp=elements[i]->U;
            std::cout << elements[i]->interpolated_u(0.0) <<  std::endl;
        }
    }
    */
    
    /* Burgers with sin
    {
        int N=100;
        double dt=0.001;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2*M_PI;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new BurgersElement();
        }
        sinconstructor(elements,minbound,maxbound);
        for (int i=0; i<2000; i++)
        {
            Galerkin_Solver_RK4_minmod(elements, dt, minbound, maxbound, 2.0);
        }
        for (int i=0; i<elements.size(); i++)
        {
            elements[i]->Utemp=elements[i]->U;
            std::cout << elements[i]->interpolated_u(0.0) << std::endl;
        }
        for (auto el : elements)
        {
            delete el;
        }
        elements.clear();
    }
    */

    /* burgers with square and minmod
    {
        int N=100;
        double dt=0.001;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new BurgersElement();
        }
        squareconstructor(elements,minbound,maxbound);
        for (int i=0; i<2000; i++)
        {
            Galerkin_Solver_RK4_minmod(elements, dt, minbound, maxbound, 2.0);
        }
        for (int i=0; i<elements.size(); i++)
        {
            elements[i]->Utemp=elements[i]->U;
            std::cout << elements[i]->interpolated_u(0.0) <<  std::endl;
        }
    }
    
    */

    /* oscillations illustration
    {
        int N=11;
        double dt=0.04;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        squareconstructor(elements,minbound,maxbound);
        for (int i=0; i<3; i++)
        {
            Galerkin_Solver(elements, dt, minbound, maxbound);
        }
        for (int i=3; i<elements.size()-2; i++)
        {
            elements[i]->Utemp=elements[i]->U;
            std::cout <<elements[i]->Utemp[0]<<  std::endl;
            std::cout <<elements[i]->Utemp[1]<<  std::endl;
        }
    }
    */
    
    // limiters illustration (rember to set square constructor back to <= and galerkin_minmod to not print)
    {
        int N=11;
        double dt=0.04;
        double t=0.0;
        double minbound=0.0;
        double maxbound=2.0;
        std::vector<AdvectionElement*> elements(N);
        for (int i = 0; i < N; ++i)
        {
            elements[i] = new AdvectionElement();
        }
        squareconstructor(elements,minbound,maxbound);
        for (int i=0; i<4; i++)
        {
            Galerkin_Solver_minmod(elements, dt, minbound, maxbound,1.0);
        }
    }

    return 0;
}