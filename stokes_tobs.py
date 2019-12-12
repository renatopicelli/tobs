from oct2py import octave
from dolfin import *
from dolfin_adjoint import *
import numpy as np

import cplex
from cplex.exceptions import CplexError

set_log_level(50)

class Optimizer(object):

    objfun_rf = None        #Objective Functions
    iter_fobj = 0
    iter_dobj = 0
    file_out = None         #Results file output
    xi_array = None         #Design Variables
    vf_fun = None
    cst_U = []
    cst_L = []
    cst_num = 0
    rho = None

    def __init__(self, fc_xi):
        self.fc_xi = Function( fc_xi.function_space(), name="Control" )
        self.nvars = len(self.fc_xi.vector())
        self.control = Control(fc_xi)

    def __check_ds_vars__(self):
        chk_var = False
        if self.xi_array is None:
            self.xi_array = np.copy(self.rho)
            chk_var = True
        else:
            xi_eval = self.xi_array - self.rho
            xi_nrm  = np.linalg.norm(xi_eval)
            if xi_nrm > 1e-16:
                self.xi_array = np.copy(self.rho)
                chk_var = True      #A variavel de projeto ja foi carregada

        if chk_var is True:
            self.fc_xi.vector()[:] = self.rho
        else:
            pass
        ds_vars = self.fc_xi
        return ds_vars

    def __vf_fun_var_assem__(self):
        fc_xi_tst   = TestFunction(self.fc_xi.function_space())
        self.vol_xi  = assemble(fc_xi_tst * Constant(1.0) * dx)
        self.vol_sum = self.vol_xi.sum()

    def add_plot_res(self, file_out):
        self.file_out = file_out

    def add_objfun(self, AD_Obj_fx):
        self.objfun_rf = ReducedFunctional(AD_Obj_fx, self.control)

    def obj_fun(self, user_data=None):
        ds_vars = self.__check_ds_vars__()
        fval = self.objfun_rf(ds_vars)
        print(" fval: ", fval)
        self.iter_fobj += 1

        return fval

    def obj_dfun(self, user_data=None):
        ds_vars  = self.__check_ds_vars__()
        self.objfun_rf(ds_vars)
        #Derivada da funcao objetivo
        dfval = self.objfun_rf.derivative().vector()
        #salva os arquivos de resultados
        if self.file_out is not None:
            self.file_out << self.fc_xi
        # contador do numero de iteracoes
        self.iter_dobj += 1

        return dfval

    def add_volf_constraint(self, upp, lwr):
        self.__vf_fun_var_assem__()

        self.cst_U.append(upp)
        self.cst_L.append(lwr)

        self.cst_num += 1

    def volfrac_fun(self):

        self.__check_ds_vars__()

        volume_val = float( self.vol_xi.inner( self.fc_xi.vector() ) )
        print('---------')
        print(volume_val/self.vol_sum)
        print('---------')

        return volume_val/self.vol_sum

    def volfrac_dfun(self, user_data=None):
        v_df = self.vol_xi/self.vol_sum
        print('*/-*/-*/-*/')
        print(np.array(v_df)[0])
        print('*/-*/-*/-*/')
        return v_df

    def flag_jacobian(self):
        rows = []
        for i in range(self.cst_num):
            rows += [i] * self.nvars
        cols = range(self.nvars) * self.cst_num

        return (np.array(rows, dtype=np.int), np.array(cols, dtype=np.int))

    def cst_fval(self, user_data=None):
        cst_val = np.array(self.volfrac_fun(), dtype=np.float)

        return cst_val.T

    def jacobian(self, flag=False, user_data=None):
        if flag:
            dfval = self.flag_jacobian()
        else:
            dfval = self.volfrac_dfun()

        return dfval

octave.addpath('~')
parameters["std_out_all_processes"] = False
pasta = "output/"

mu = Constant(1.0)                   # viscosity
alphaunderbar = 2.5 * mu / (100**2)  # parameter for \alpha
alphabar = 2.5 * mu / (0.01**2)      # parameter for \alpha
alphaunderbar = 2.5 * mu *1.e-10
alphabar = 2.5 * mu * 1e10
q = Constant(1.0) # q value that controls difficulty/discrete-valuedness of solution

def alpha(rho):
    """Inverse permeability as a function of rho, equation (40)"""
    return alphabar + (alphaunderbar - alphabar) * rho * (1 + q) / (rho + q)

N = 40
delta = 1.5  # The aspect ratio of the domain, 1 high and \delta wide
V = Constant(1.0/3) * delta  # want the fluid to occupy 1/3 of the domain

mesh = Mesh(RectangleMesh(Point(0.0, 0.0), Point(delta, 1.0), int(N*3/2), N, diagonal="crossed"))
A = FunctionSpace(mesh, "DG", 0)        # control function space

U_h = VectorElement("CG", mesh.ufl_cell(), 2)
P_h = FiniteElement("CG", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, U_h*P_h)          # mixed Taylor-Hood function space

# Define the boundary condition on velocity

class InflowOutflow(UserExpression):
    def eval(self, values, x):
        values[1] = 0.0
        values[0] = 0.0
        l = 1.0/6.0
        gbar = 1.0

        if x[0] == 0.0 or x[0] == delta:
            if (1.0/4 - l/2) < x[1] < (1.0/4 + l/2):
                t = x[1] - 1.0/4
                values[0] = gbar*(1 - (2*t/l)**2)
            if (3.0/4 - l/2) < x[1] < (3.0/4 + l/2):
                t = x[1] - 3.0/4
                values[0] = gbar*(1 - (2*t/l)**2)

    def value_shape(self):
        return (2,)

def forward(rho):
    """Solve the forward problem for a given fluid distribution rho(x)."""
    w_resp = Function(W)
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)

    F = (alpha(rho) * inner(u, v) * dx + mu*inner(grad(u)+grad(u).T, grad(v)) * dx +
         inner(grad(p), v) * dx  + inner(div(u), q) * dx)
    bc = DirichletBC(W.sub(0), InflowOutflow(degree=1), "on_boundary")
    solve(lhs(F) == rhs(F), w_resp, bcs=bc)

    return w_resp

class Distribution(UserExpression):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def eval_cell(self, values, x, ufc_cell):
        values[0] = 1

    def value_shape(self):
        return ()

def cplex_optimize(prob, nvar, my_obj, my_constcoef, my_rlimits, my_ll, my_ul):
    prob.objective.set_sense(prob.objective.sense.minimize)

    my_ctype = "I"*nvar
    my_colnames = ["x"+str(item) for item in range(nvar)]
    my_sense = ["L", "G"]
    my_rownames = ["r1", "r2"]

    prob.variables.add(obj=my_obj, lb=my_ll, ub=my_ul, types=my_ctype,
                       names=my_colnames)

    rows = [cplex.SparsePair(ind=["x"+str(item) for item in range(nvar)], val = my_constcoef[0]),
            cplex.SparsePair(ind=["x"+str(item) for item in range(nvar)], val = my_constcoef[1])]

    prob.linear_constraints.add(lin_expr=rows, senses=my_sense, rhs=my_rlimits, names=my_rownames)


if __name__ == "__main__":
    rho = interpolate(Distribution(), A)
    rho.rename("control", "")
    w_resp   = forward(rho)
    (u, p) = w_resp.split()

    controls = File(pasta + "control.pvd")
    state_file = File(pasta + "veloc.pvd")
    rho_viz = Function(A, name="ControlVisualisation")

    controls << rho

    iteration = 0
    epsilons = .2

    while True:
        J = assemble(0.5 * inner(alpha(rho) * u, u) * dx + mu * inner(grad(u), grad(u)) * dx)
        nvar = len(rho.vector())

        fval = Optimizer(rho)
        fval.add_objfun(J)
        fval.add_volf_constraint(0.7,0.5)
        fval.rho = rho.vector()
        fid_rho = File("rho.pvd")
        x_L = np.ones((nvar), dtype=np.float) * 0.0
        x_U = np.ones((nvar), dtype=np.float) * 1.0
        acst_L = np.array(fval.cst_L)
        acst_U = np.array(fval.cst_U)
        j = float(fval.obj_fun(rho.vector()))
        if iteration == 0: jd_previous = np.array(fval.obj_dfun()).reshape((-1,1))
        jd = (np.array(fval.obj_dfun()).reshape((-1,1)) + jd_previous)/2 #stabilization
        cs = fval.cst_fval()
        jac = np.array(fval.jacobian()).reshape((-1,1))
        ans = octave.stokes(
                nvar,
                x_L,
                x_U,
                fval.cst_num,
                acst_L,
                acst_U,
                j,
                jd,
                cs,
                jac,
                iteration,
                epsilons,
                np.array(rho.vector())
                )
        PythonObjCoeff = ans[0][1] #because [0][0] is the design variable
        PythonConstCoeff = ans[0][2]
        PythonRelaxedLimits = ans[0][3]
        PythonLowerLimits = ans[0][4]
        PythonUpperLimits = ans[0][5]
        PythonnDesignVariables = ans[0][6]
        my_prob = cplex.Cplex()
        coef = [item[0] for item in PythonObjCoeff.tolist()]
        constcoef = PythonConstCoeff.tolist()
        rlimits = [item[0] for item in PythonRelaxedLimits.tolist()]
        ll = [item[0] for item in PythonLowerLimits.tolist()]
        ul = [item[0] for item in PythonUpperLimits.tolist()]
        cplex_optimize(my_prob, nvar, coef, constcoef, rlimits, ll, ul)

        my_prob.solve()
        design_variables = my_prob.solution.get_values()

        rho.rename("control", "")
        rho.vector().add_local(np.array(design_variables))
        controls << rho

        w_resp   = forward(rho)
        (u, p) = w_resp.split()
        u.rename("Velocidade", "")
        state_file << u
        if iteration == 700:
            break
        else: iteration += 1
        jd_previous = jd
        del my_prob

