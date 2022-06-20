// https://observablehq.com/@kelleyvanevert/katex-within-markdown@404
import define1 from "./f0bfa03e8410e4c7@404.js";

function _1(md){return(
md`# KaTeX within markdown

By default, if you want LaTeX intermingled within Markdown, you have to write:

~~~
md\`mdtext \${tex\`latext\`} mdtext\`
~~~

If you do that a lot, you may prefer the following syntax:

~~~
texmd\`mdtext \$latext\$ mdtext\`
~~~

This work was started by [Kelley van Evert](https://observablehq.com/@kelleyvanevert/katex-within-markdown). That implementation has one major problem, which is that \`\${}\` interpolation inside \`texmd\` can't return a DOM element:
`
)}

function _2(html,texmd_upstream)
{
  function makeelt(str) {
    return html`<span>${str}</span>`
  }

  return texmd_upstream`$\beta$ This *should* be useful: ${makeelt("this *isn't* markdown")} and it ${makeelt("*isn't* $te_x$")}, ${"boo"}!`
}


function _3(md){return(
md`This is more like what you might expect:`
)}

function _4(html,texmd)
{
  function makeelt(str) {
    return html`<span>${str}</span>`
  }

  return texmd`$\beta$ This is *more* useful: ${makeelt("this *isn't* markdown")} and it ${makeelt("*isn't* $te_x$")}, ${"yay"}!`
}


function _5(md){return(
md`## Example`
)}

function _6(texmd){return(
texmd`Let $\Gamma$ be a $D$-set, $A : \Gamma \to D$**-set**, and $B:\sigma(\Gamma,A) \to D$**-set**. Define the $D$-set element for internal relativized evaluation of an internal dependent product:

$$ \begin{aligned} {\sf{eval}}_{\Gamma,A,B} &\in \pi\big(\,\pi_\Gamma(A,B) \times \pi(\Gamma,A),\, n \mapsto \pi(\Gamma,B(\langle -, n_1(-) \rangle)) \,\big) \\ {\sf{eval}}_{\Gamma,A,B}&(n)(\gamma) \triangleq n_0(\gamma)(n_1(\gamma)) \end{aligned} $$

This definition is valid, because $ {\sf{eval}}_{\Gamma,A,B}(p,a) \in \pi(\Gamma,B(\langle-,a(-)\rangle)) $ is tracked by $\Lambda z.\,(p \cdot z)\cdot(a \cdot z)$, and subsequently $ {\sf{eval}}_{\Gamma,A,B} $ is tracked by $\Lambda n.\Lambda z.\,(n_0 \cdot z)\cdot(n_1 \cdot z)$.

Informally, this is of course just a messy internal $D$**-set** equivalent of the following picture: $$ {\sf{eval}}_{\Gamma,A,B} : {\big({B_{(\gamma,a)}}^{A_\gamma}\big)}^{\Gamma} \times {A_\gamma}^\Gamma \to {B_{(\gamma,a)}}^\Gamma$$`
)}

function _7(md){return(
md`## Usage`
)}

function _8(md){return(
md`Use <code>texmd</code> as you would <code>tex</code>, <code>md</code> or <code>html</code>.

You can use inline maths $ ... $ as well as display maths $$ ... $$.`
)}

function _9(md){return(
md`## Implementation`
)}

function _texmd(md,tex){return(
function texmd() {
  // First perform tex replacement on all the strings in the template
  // Each piece of tex is replaced with a placeholder for now, to keep it away from the markdown parser
  var display_maths = [];
  var inline_maths = [];
  
  // Note use of raw strings here; important to keep backslashes intact so tex can see them
  let strings = arguments[0].raw.map(function(string) {
    // Then extract all the TeX code and replace it with placeholders (so it won't be run through md)
    string = string.replace(/\$\$([^\$]*)\$\$/g, (s, m) => {
      var i = display_maths.push(m) - 1;
      return `\n\n<div class="dmath" n="${i}"></div>\n\n`;
    });

    string = string.replace(/\$([^\$]*)\$/g, (s, m) => {
      var i = inline_maths.push(m) - 1;
      return `<span class="imath" n="${i}"></span>`;
    });

    return string
  })
  
  // Now we can run this through markdown with proper ${} substitutions
  var node = md(strings, ...Array.prototype.slice.call(arguments, 1));
  
  // And finally replace TeX placeholders with output of tex()
  node.querySelectorAll('span.imath').forEach(span => {
    var i = parseInt(span.attributes.n.value);
    if (inline_maths[i] === undefined) throw 'What?';
    span.appendChild(tex`${inline_maths[i]}`);
    span.style['font-size'] = '.95em';
  });
  
  node.querySelectorAll('div.dmath').forEach(div => {
    var i = parseInt(div.attributes.n.value);
    div.appendChild(tex.block`${display_maths[i]}`);
    div.style['font-size'] = '.95em';
  });
  
  return node;
}
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], _1);
  main.variable(observer()).define(["html","texmd_upstream"], _2);
  main.variable(observer()).define(["md"], _3);
  main.variable(observer()).define(["html","texmd"], _4);
  main.variable(observer()).define(["md"], _5);
  main.variable(observer()).define(["texmd"], _6);
  main.variable(observer()).define(["md"], _7);
  main.variable(observer()).define(["md"], _8);
  main.variable(observer()).define(["md"], _9);
  main.variable(observer("texmd")).define("texmd", ["md","tex"], _texmd);
  const child1 = runtime.module(define1);
  main.import("texmd", "texmd_upstream", child1);
  return main;
}
