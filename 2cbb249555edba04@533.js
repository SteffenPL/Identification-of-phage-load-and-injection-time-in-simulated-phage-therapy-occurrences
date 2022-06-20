// https://observablehq.com/@3784219e03ed337e/interactive-phage-simulation@533
import define1 from "./f0bfa03e8410e4c7@404.js";

function _1(md){return(
md`# Identification of phage load and injection time in simulated phage therapy occurrences
## by Steffen Plunder and Luigi Marongiu
Available at: https://observablehq.com/@3784219e03ed337e/interactive-phage-simulation

### Supplementary data: dynamic plot
Dynamic code for case 1 reported in the paper. The scenario refers to the data reported by Payne, R.J. and Jansen, V.A. "Understanding Bacteriophage Therapy as a Density-Dependent Kinetic Process", *J. Theor. Biol. 2001, 208, 37–48*.

The default values for the parameters below as well as the ODEs are set according to the paper. The output of this computation is the plot of the densities of the hypothetical bacterial and phagial populations. The plot has been extended five-fold with respect to the work of Payne and Jansen (100 hours versus 20 hours). 

The user can have a feeling on how phage therapy depends on case-specific conditions by modifying the parameters using the slider below. The list of parameters is as follows:

* tᵢₙ = time of administration of the phage dose   (h)
* vᵢₙ = amount of the phage dose                   (PFU/mL)
* S₀  = initial density of susceptible bacteria    (CFU/mL)
* K = carrying capacity of the environment         (CFU/mL)
* μ = maximum growth rate of the bacteria          (1/h)
* ω = outflow of the environment                   (1/mL)
* δ = infection rate of the phage                  (1/min)
* η = inverse of latency time                      (min)
* λ = degradation rate of the phage                (1/day)
* β = burst size of the phage                      (PFU)`
)}

function _T_in(Inputs,texmd){return(
Inputs.range([0, 60], {value: 15.4, step: 0.1, label: texmd`$t_{in}$`})
)}

function _log_v_in(Inputs,texmd){return(
Inputs.range([0, 10], 
                           {value: 8.8, step: 0.1,  
                            label: texmd`$\log_{10}(v_{in})$`})
)}

function _log_su0(Inputs,texmd){return(
Inputs.range([0, 9], {value: 3, step: 0.1, label: texmd`$\log_{10}(S_{0})$`})
)}

function _log_K(Inputs,texmd){return(
Inputs.range([0, 9], {value: 7, step: 0.1, label: texmd`$\log_{10}(K)$`})
)}

function _Mu(Inputs,texmd){return(
Inputs.range([0, 1], {value: 0.5, step: 0.1, label: texmd`$\mu$`})
)}

function _Omega(Inputs,texmd){return(
Inputs.range([0, 1], {value: 0.2, step: 0.1, label: texmd`$\omega$`})
)}

function _char_delta(Inputs,texmd){return(
Inputs.range([0.1, 10], {value: 1.0, step: 0.1, label: texmd`$\delta$`})
)}

function _Eta(Inputs,texmd){return(
Inputs.range([0, 10], {value: 5, step: 0.1, label: texmd`$\eta$`})
)}

function _Lambda(Inputs,texmd){return(
Inputs.range([0, 10], {value: 5, step: 0.1, label: texmd`$\lambda$`})
)}

function _Beta(Inputs,texmd){return(
Inputs.range([1, 1000], {value: 100, step: 0.1, label: texmd`$\beta$`})
)}

function _12(Plot,names){return(
Plot.legend({color: {type: "ordinal", domain: names, range: ["green", "red", "blue"]}})
)}

function _13(Plot,case1,sol,terminal_time,p){return(
Plot.plot({
  y: {type: "log", domain: [1e2, 1e8], label: "↑ Concentration"},
  x: {domain: [0.0, case1.t_end], label: "Time (in h)→"},
  width: 1000,
  height: 400,
  grid: true,
  color: {
    legend: true
  },
  marks: [
    Plot.line(sol, {x: "t", y: "u1", stroke: "green"}),    
    Plot.line(sol, {x: "t", y: "u2", stroke: "red"}),    
    Plot.line(sol, {x: "t", y: "u3", stroke: "blue"}),
    Plot.link(terminal_time, 
      {
        x1:"t", 
        x2:"t",  
        y1:"min", 
        y2:"max", 
        stroke: "gray", 
        strokeWidth: 1
      }),    
    Plot.text(terminal_time, 
      {
      x: "t",
      y: "text_y",
      text: "text",
      fill: "currentColor",
      stroke: "white",
      dx: 40
    }),
    Plot.link([p], 
      {
        x1:"t_in", 
        x2:"t_in",  
        y1: "t_in", 
        y2:"v_in", 
        stroke: "blue", 
        strokeWidth: 1
      }),    
  ]
})
)}

function _14(md){return(
md`### Computations:`
)}

function _extinction_condition(Inputs){return(
Inputs.select(["yes", "no"], {label: "Extinction condition:"})
)}

function _dt(){return(
0.01
)}

function _Kappa(log_K){return(
Math.pow(10, log_K)
)}

function _Delta(char_delta){return(
char_delta*Math.pow(10, -7)
)}

function _S0(log_su0){return(
Math.pow(10, log_su0)
)}

function _V_in(log_v_in){return(
Math.pow(10, log_v_in)
)}

function _names(){return(
["uninfected bacteria", "lytic bacteria", "free phage"]
)}

function _terminal_time(){return(
[{t: 20, min:1, max:1e10, text: "terminal time", text_y: 1e8}]
)}

function _f_logistic(){return(
function f_logistic(x, p)
{ 
  //const N = x.u1 + x.u2 + x.u3
  const N = x.u1 + x.u2 
  const rho = 1 - N/p.kappa
  
  return {
    //u1: (p.mu*x.u1*rho) - (p.delta*x.u1*x.u4)                    - p.omega*x.u1,
    //u2: (p.mu*x.u2*rho) + (p.delta*x.u1*x.u4) - (p.eta*x.u2)     - p.omega*x.u2,
    //u3: (p.nu*x.u3*rho)                                          - p.omega*x.u3,
    //u4: (p.mu*p.beta*x.u2) - (p.delta*x.u1*x.u4) - (p.lambda*x.u4) - p.omega*x.u1

    u1: (p.mu*x.u1*rho)     - (p.delta*x.u1*x.u3)                   - p.omega*x.u1,
    u2: (p.mu*x.u2*rho)     + (p.delta*x.u1*x.u3) - (p.eta*x.u2)    - p.omega*x.u2,
    u3: (p.eta*p.beta*x.u2) - (p.delta*x.u1*x.u3) - (p.lambda*x.u3) - p.omega*x.u3
  }
}
)}

function _compute_solution(extinction_condition,f_logistic){return(
function compute_solution(ode, x0, t_end, p, dt = 0.01)
{
  var t = 0.0 
  const n_steps = Math.ceil(t_end/dt)

  const extinc = (extinction_condition == "yes")
  
  var x = Object.assign({}, x0)  // copy object
  x.t = t
  var sol = []
  sol[0] = Object.assign({}, x0)

  for (let k = 1; k < n_steps; ++k) {
    x.t += dt

    // handle the special event of phage injection
    if( x.t + dt >= p.t_in && x.t < p.t_in ){
      x.u3 += p.v_in
    }

    var du = f_logistic(x, p)


    // handel is one population is below 1
    if( extinc )
    {
      if(  x.u1 < 1.0 && du.u1 < 0 ){
        x.u1 = 0.0
      }
      
      if(  x.u2 < 1.0 && du.u2 < 0 ){
        x.u2 = 0.0
      }
      
      if(  x.u3 < 1.0 && du.u3 < 0 ){
        x.u3 = 0.0
      }
    }
    
    
    x.u1 += dt * du.u1
    x.u2 += dt * du.u2
    x.u3 += dt * du.u3
   // x.u4 += dt * du.u4
    
    sol[k] = Object.assign({}, x)
  }

  return sol
}
)}

function _case1(S0,V_in,T_in,Omega,Kappa,Mu,Delta,Eta,Lambda,Beta)
{
    const su0   = S0       // initial susceptible population
    const v0    = 0        // initial phage population
    const v_in  = V_in     // phage inoculum
    const t_in  = T_in     // time of phage inoculum
    const t_mx  = 100.0    // time span 0-tmax
    const omega = Omega    // outflow (/h)
    const kappa = Kappa    // carrying capacity (based on M, case 1)
    const mu     = Mu      // maximum growth rate susceptible (/h)
    const delta = Delta    // adsorption rate       [-]
    const eta   = Eta      // reciprocal of latency [1/h]
    const lambda = Lambda  // decay rate in hours   [PFU/h]
    const beta  = Beta     // burst size            [PFU]

    return {x0: {u1:su0, u2:0.0, u3:v0}, p: {mu, lambda, delta, beta, eta, kappa, omega, t_in, v_in}, t_end: t_mx }
}


function _p(case1,T_in,V_in)
{
  // add the parameters from the slides to the case. This is done with the function Object.assign which merges two objectss.
  const p = case1.p
  p.t_in = T_in 
  p.v_in = Math.pow(10, V_in)

  return Object.assign({}, case1.p, {v_in: V_in, t_in: T_in})
}


function _sol(compute_solution,f_logistic,case1,p,dt){return(
compute_solution(f_logistic, case1.x0, case1.t_end, p, dt)
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], _1);
  main.variable(observer("viewof T_in")).define("viewof T_in", ["Inputs","texmd"], _T_in);
  main.variable(observer("T_in")).define("T_in", ["Generators", "viewof T_in"], (G, _) => G.input(_));
  main.variable(observer("viewof log_v_in")).define("viewof log_v_in", ["Inputs","texmd"], _log_v_in);
  main.variable(observer("log_v_in")).define("log_v_in", ["Generators", "viewof log_v_in"], (G, _) => G.input(_));
  main.variable(observer("viewof log_su0")).define("viewof log_su0", ["Inputs","texmd"], _log_su0);
  main.variable(observer("log_su0")).define("log_su0", ["Generators", "viewof log_su0"], (G, _) => G.input(_));
  main.variable(observer("viewof log_K")).define("viewof log_K", ["Inputs","texmd"], _log_K);
  main.variable(observer("log_K")).define("log_K", ["Generators", "viewof log_K"], (G, _) => G.input(_));
  main.variable(observer("viewof Mu")).define("viewof Mu", ["Inputs","texmd"], _Mu);
  main.variable(observer("Mu")).define("Mu", ["Generators", "viewof Mu"], (G, _) => G.input(_));
  main.variable(observer("viewof Omega")).define("viewof Omega", ["Inputs","texmd"], _Omega);
  main.variable(observer("Omega")).define("Omega", ["Generators", "viewof Omega"], (G, _) => G.input(_));
  main.variable(observer("viewof char_delta")).define("viewof char_delta", ["Inputs","texmd"], _char_delta);
  main.variable(observer("char_delta")).define("char_delta", ["Generators", "viewof char_delta"], (G, _) => G.input(_));
  main.variable(observer("viewof Eta")).define("viewof Eta", ["Inputs","texmd"], _Eta);
  main.variable(observer("Eta")).define("Eta", ["Generators", "viewof Eta"], (G, _) => G.input(_));
  main.variable(observer("viewof Lambda")).define("viewof Lambda", ["Inputs","texmd"], _Lambda);
  main.variable(observer("Lambda")).define("Lambda", ["Generators", "viewof Lambda"], (G, _) => G.input(_));
  main.variable(observer("viewof Beta")).define("viewof Beta", ["Inputs","texmd"], _Beta);
  main.variable(observer("Beta")).define("Beta", ["Generators", "viewof Beta"], (G, _) => G.input(_));
  main.variable(observer()).define(["Plot","names"], _12);
  main.variable(observer()).define(["Plot","case1","sol","terminal_time","p"], _13);
  main.variable(observer()).define(["md"], _14);
  main.variable(observer("viewof extinction_condition")).define("viewof extinction_condition", ["Inputs"], _extinction_condition);
  main.variable(observer("extinction_condition")).define("extinction_condition", ["Generators", "viewof extinction_condition"], (G, _) => G.input(_));
  main.variable(observer("dt")).define("dt", _dt);
  main.variable(observer("Kappa")).define("Kappa", ["log_K"], _Kappa);
  main.variable(observer("Delta")).define("Delta", ["char_delta"], _Delta);
  main.variable(observer("S0")).define("S0", ["log_su0"], _S0);
  main.variable(observer("V_in")).define("V_in", ["log_v_in"], _V_in);
  main.variable(observer("names")).define("names", _names);
  main.variable(observer("terminal_time")).define("terminal_time", _terminal_time);
  main.variable(observer("f_logistic")).define("f_logistic", _f_logistic);
  main.variable(observer("compute_solution")).define("compute_solution", ["extinction_condition","f_logistic"], _compute_solution);
  main.variable(observer("case1")).define("case1", ["S0","V_in","T_in","Omega","Kappa","Mu","Delta","Eta","Lambda","Beta"], _case1);
  main.variable(observer("p")).define("p", ["case1","T_in","V_in"], _p);
  main.variable(observer("sol")).define("sol", ["compute_solution","f_logistic","case1","p","dt"], _sol);
  const child1 = runtime.module(define1);
  main.import("texmd", child1);
  return main;
}
