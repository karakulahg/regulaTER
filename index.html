<!DOCTYPE html>

<html>

<head>

<meta charset="utf-8" />
<meta name="generator" content="pandoc" />
<meta http-equiv="X-UA-Compatible" content="IE=EDGE" />

<meta name="viewport" content="width=device-width, initial-scale=1" />



<title>regulaTER: an R package for analysis of transposable elements in accessible regions</title>

<script src="data:application/javascript;base64,Ly8gUGFuZG9jIDIuOSBhZGRzIGF0dHJpYnV0ZXMgb24gYm90aCBoZWFkZXIgYW5kIGRpdi4gV2UgcmVtb3ZlIHRoZSBmb3JtZXIgKHRvCi8vIGJlIGNvbXBhdGlibGUgd2l0aCB0aGUgYmVoYXZpb3Igb2YgUGFuZG9jIDwgMi44KS4KZG9jdW1lbnQuYWRkRXZlbnRMaXN0ZW5lcignRE9NQ29udGVudExvYWRlZCcsIGZ1bmN0aW9uKGUpIHsKICB2YXIgaHMgPSBkb2N1bWVudC5xdWVyeVNlbGVjdG9yQWxsKCJkaXYuc2VjdGlvbltjbGFzcyo9J2xldmVsJ10gPiA6Zmlyc3QtY2hpbGQiKTsKICB2YXIgaSwgaCwgYTsKICBmb3IgKGkgPSAwOyBpIDwgaHMubGVuZ3RoOyBpKyspIHsKICAgIGggPSBoc1tpXTsKICAgIGlmICghL15oWzEtNl0kL2kudGVzdChoLnRhZ05hbWUpKSBjb250aW51ZTsgIC8vIGl0IHNob3VsZCBiZSBhIGhlYWRlciBoMS1oNgogICAgYSA9IGguYXR0cmlidXRlczsKICAgIHdoaWxlIChhLmxlbmd0aCA+IDApIGgucmVtb3ZlQXR0cmlidXRlKGFbMF0ubmFtZSk7CiAgfQp9KTsK"></script>

<style type="text/css">
  code{white-space: pre-wrap;}
  span.smallcaps{font-variant: small-caps;}
  span.underline{text-decoration: underline;}
  div.column{display: inline-block; vertical-align: top; width: 50%;}
  div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
  ul.task-list{list-style: none;}
    </style>


<style type="text/css">
  code {
    white-space: pre;
  }
  .sourceCode {
    overflow: visible;
  }
</style>
<style type="text/css" data-origin="pandoc">
pre > code.sourceCode { white-space: pre; position: relative; }
pre > code.sourceCode > span { display: inline-block; line-height: 1.25; }
pre > code.sourceCode > span:empty { height: 1.2em; }
.sourceCode { overflow: visible; }
code.sourceCode > span { color: inherit; text-decoration: inherit; }
div.sourceCode { margin: 1em 0; }
pre.sourceCode { margin: 0; }
@media screen {
div.sourceCode { overflow: auto; }
}
@media print {
pre > code.sourceCode { white-space: pre-wrap; }
pre > code.sourceCode > span { text-indent: -5em; padding-left: 5em; }
}
pre.numberSource code
  { counter-reset: source-line 0; }
pre.numberSource code > span
  { position: relative; left: -4em; counter-increment: source-line; }
pre.numberSource code > span > a:first-child::before
  { content: counter(source-line);
    position: relative; left: -1em; text-align: right; vertical-align: baseline;
    border: none; display: inline-block;
    -webkit-touch-callout: none; -webkit-user-select: none;
    -khtml-user-select: none; -moz-user-select: none;
    -ms-user-select: none; user-select: none;
    padding: 0 4px; width: 4em;
    color: #aaaaaa;
  }
pre.numberSource { margin-left: 3em; border-left: 1px solid #aaaaaa;  padding-left: 4px; }
div.sourceCode
  {   }
@media screen {
pre > code.sourceCode > span > a:first-child::before { text-decoration: underline; }
}
code span.al { color: #ff0000; font-weight: bold; } /* Alert */
code span.an { color: #60a0b0; font-weight: bold; font-style: italic; } /* Annotation */
code span.at { color: #7d9029; } /* Attribute */
code span.bn { color: #40a070; } /* BaseN */
code span.bu { } /* BuiltIn */
code span.cf { color: #007020; font-weight: bold; } /* ControlFlow */
code span.ch { color: #4070a0; } /* Char */
code span.cn { color: #880000; } /* Constant */
code span.co { color: #60a0b0; font-style: italic; } /* Comment */
code span.cv { color: #60a0b0; font-weight: bold; font-style: italic; } /* CommentVar */
code span.do { color: #ba2121; font-style: italic; } /* Documentation */
code span.dt { color: #902000; } /* DataType */
code span.dv { color: #40a070; } /* DecVal */
code span.er { color: #ff0000; font-weight: bold; } /* Error */
code span.ex { } /* Extension */
code span.fl { color: #40a070; } /* Float */
code span.fu { color: #06287e; } /* Function */
code span.im { } /* Import */
code span.in { color: #60a0b0; font-weight: bold; font-style: italic; } /* Information */
code span.kw { color: #007020; font-weight: bold; } /* Keyword */
code span.op { color: #666666; } /* Operator */
code span.ot { color: #007020; } /* Other */
code span.pp { color: #bc7a00; } /* Preprocessor */
code span.sc { color: #4070a0; } /* SpecialChar */
code span.ss { color: #bb6688; } /* SpecialString */
code span.st { color: #4070a0; } /* String */
code span.va { color: #19177c; } /* Variable */
code span.vs { color: #4070a0; } /* VerbatimString */
code span.wa { color: #60a0b0; font-weight: bold; font-style: italic; } /* Warning */

</style>
<script>
// apply pandoc div.sourceCode style to pre.sourceCode instead
(function() {
  var sheets = document.styleSheets;
  for (var i = 0; i < sheets.length; i++) {
    if (sheets[i].ownerNode.dataset["origin"] !== "pandoc") continue;
    try { var rules = sheets[i].cssRules; } catch (e) { continue; }
    for (var j = 0; j < rules.length; j++) {
      var rule = rules[j];
      // check if there is a div.sourceCode rule
      if (rule.type !== rule.STYLE_RULE || rule.selectorText !== "div.sourceCode") continue;
      var style = rule.style.cssText;
      // check if color or background-color is set
      if (rule.style.color === '' && rule.style.backgroundColor === '') continue;
      // replace div.sourceCode by a pre.sourceCode rule
      sheets[i].deleteRule(j);
      sheets[i].insertRule('pre.sourceCode{' + style + '}', j);
    }
  }
})();
</script>




<link rel="stylesheet" href="data:text/css,body%20%7B%0Abackground%2Dcolor%3A%20%23fff%3B%0Amargin%3A%201em%20auto%3B%0Amax%2Dwidth%3A%20700px%3B%0Aoverflow%3A%20visible%3B%0Apadding%2Dleft%3A%202em%3B%0Apadding%2Dright%3A%202em%3B%0Afont%2Dfamily%3A%20%22Open%20Sans%22%2C%20%22Helvetica%20Neue%22%2C%20Helvetica%2C%20Arial%2C%20sans%2Dserif%3B%0Afont%2Dsize%3A%2014px%3B%0Aline%2Dheight%3A%201%2E35%3B%0A%7D%0A%23TOC%20%7B%0Aclear%3A%20both%3B%0Amargin%3A%200%200%2010px%2010px%3B%0Apadding%3A%204px%3B%0Awidth%3A%20400px%3B%0Aborder%3A%201px%20solid%20%23CCCCCC%3B%0Aborder%2Dradius%3A%205px%3B%0Abackground%2Dcolor%3A%20%23f6f6f6%3B%0Afont%2Dsize%3A%2013px%3B%0Aline%2Dheight%3A%201%2E3%3B%0A%7D%0A%23TOC%20%2Etoctitle%20%7B%0Afont%2Dweight%3A%20bold%3B%0Afont%2Dsize%3A%2015px%3B%0Amargin%2Dleft%3A%205px%3B%0A%7D%0A%23TOC%20ul%20%7B%0Apadding%2Dleft%3A%2040px%3B%0Amargin%2Dleft%3A%20%2D1%2E5em%3B%0Amargin%2Dtop%3A%205px%3B%0Amargin%2Dbottom%3A%205px%3B%0A%7D%0A%23TOC%20ul%20ul%20%7B%0Amargin%2Dleft%3A%20%2D2em%3B%0A%7D%0A%23TOC%20li%20%7B%0Aline%2Dheight%3A%2016px%3B%0A%7D%0Atable%20%7B%0Amargin%3A%201em%20auto%3B%0Aborder%2Dwidth%3A%201px%3B%0Aborder%2Dcolor%3A%20%23DDDDDD%3B%0Aborder%2Dstyle%3A%20outset%3B%0Aborder%2Dcollapse%3A%20collapse%3B%0A%7D%0Atable%20th%20%7B%0Aborder%2Dwidth%3A%202px%3B%0Apadding%3A%205px%3B%0Aborder%2Dstyle%3A%20inset%3B%0A%7D%0Atable%20td%20%7B%0Aborder%2Dwidth%3A%201px%3B%0Aborder%2Dstyle%3A%20inset%3B%0Aline%2Dheight%3A%2018px%3B%0Apadding%3A%205px%205px%3B%0A%7D%0Atable%2C%20table%20th%2C%20table%20td%20%7B%0Aborder%2Dleft%2Dstyle%3A%20none%3B%0Aborder%2Dright%2Dstyle%3A%20none%3B%0A%7D%0Atable%20thead%2C%20table%20tr%2Eeven%20%7B%0Abackground%2Dcolor%3A%20%23f7f7f7%3B%0A%7D%0Ap%20%7B%0Amargin%3A%200%2E5em%200%3B%0A%7D%0Ablockquote%20%7B%0Abackground%2Dcolor%3A%20%23f6f6f6%3B%0Apadding%3A%200%2E25em%200%2E75em%3B%0A%7D%0Ahr%20%7B%0Aborder%2Dstyle%3A%20solid%3B%0Aborder%3A%20none%3B%0Aborder%2Dtop%3A%201px%20solid%20%23777%3B%0Amargin%3A%2028px%200%3B%0A%7D%0Adl%20%7B%0Amargin%2Dleft%3A%200%3B%0A%7D%0Adl%20dd%20%7B%0Amargin%2Dbottom%3A%2013px%3B%0Amargin%2Dleft%3A%2013px%3B%0A%7D%0Adl%20dt%20%7B%0Afont%2Dweight%3A%20bold%3B%0A%7D%0Aul%20%7B%0Amargin%2Dtop%3A%200%3B%0A%7D%0Aul%20li%20%7B%0Alist%2Dstyle%3A%20circle%20outside%3B%0A%7D%0Aul%20ul%20%7B%0Amargin%2Dbottom%3A%200%3B%0A%7D%0Apre%2C%20code%20%7B%0Abackground%2Dcolor%3A%20%23f7f7f7%3B%0Aborder%2Dradius%3A%203px%3B%0Acolor%3A%20%23333%3B%0Awhite%2Dspace%3A%20pre%2Dwrap%3B%20%0A%7D%0Apre%20%7B%0Aborder%2Dradius%3A%203px%3B%0Amargin%3A%205px%200px%2010px%200px%3B%0Apadding%3A%2010px%3B%0A%7D%0Apre%3Anot%28%5Bclass%5D%29%20%7B%0Abackground%2Dcolor%3A%20%23f7f7f7%3B%0A%7D%0Acode%20%7B%0Afont%2Dfamily%3A%20Consolas%2C%20Monaco%2C%20%27Courier%20New%27%2C%20monospace%3B%0Afont%2Dsize%3A%2085%25%3B%0A%7D%0Ap%20%3E%20code%2C%20li%20%3E%20code%20%7B%0Apadding%3A%202px%200px%3B%0A%7D%0Adiv%2Efigure%20%7B%0Atext%2Dalign%3A%20center%3B%0A%7D%0Aimg%20%7B%0Abackground%2Dcolor%3A%20%23FFFFFF%3B%0Apadding%3A%202px%3B%0Aborder%3A%201px%20solid%20%23DDDDDD%3B%0Aborder%2Dradius%3A%203px%3B%0Aborder%3A%201px%20solid%20%23CCCCCC%3B%0Amargin%3A%200%205px%3B%0A%7D%0Ah1%20%7B%0Amargin%2Dtop%3A%200%3B%0Afont%2Dsize%3A%2035px%3B%0Aline%2Dheight%3A%2040px%3B%0A%7D%0Ah2%20%7B%0Aborder%2Dbottom%3A%204px%20solid%20%23f7f7f7%3B%0Apadding%2Dtop%3A%2010px%3B%0Apadding%2Dbottom%3A%202px%3B%0Afont%2Dsize%3A%20145%25%3B%0A%7D%0Ah3%20%7B%0Aborder%2Dbottom%3A%202px%20solid%20%23f7f7f7%3B%0Apadding%2Dtop%3A%2010px%3B%0Afont%2Dsize%3A%20120%25%3B%0A%7D%0Ah4%20%7B%0Aborder%2Dbottom%3A%201px%20solid%20%23f7f7f7%3B%0Amargin%2Dleft%3A%208px%3B%0Afont%2Dsize%3A%20105%25%3B%0A%7D%0Ah5%2C%20h6%20%7B%0Aborder%2Dbottom%3A%201px%20solid%20%23ccc%3B%0Afont%2Dsize%3A%20105%25%3B%0A%7D%0Aa%20%7B%0Acolor%3A%20%230033dd%3B%0Atext%2Ddecoration%3A%20none%3B%0A%7D%0Aa%3Ahover%20%7B%0Acolor%3A%20%236666ff%3B%20%7D%0Aa%3Avisited%20%7B%0Acolor%3A%20%23800080%3B%20%7D%0Aa%3Avisited%3Ahover%20%7B%0Acolor%3A%20%23BB00BB%3B%20%7D%0Aa%5Bhref%5E%3D%22http%3A%22%5D%20%7B%0Atext%2Ddecoration%3A%20underline%3B%20%7D%0Aa%5Bhref%5E%3D%22https%3A%22%5D%20%7B%0Atext%2Ddecoration%3A%20underline%3B%20%7D%0A%0Acode%20%3E%20span%2Ekw%20%7B%20color%3A%20%23555%3B%20font%2Dweight%3A%20bold%3B%20%7D%20%0Acode%20%3E%20span%2Edt%20%7B%20color%3A%20%23902000%3B%20%7D%20%0Acode%20%3E%20span%2Edv%20%7B%20color%3A%20%2340a070%3B%20%7D%20%0Acode%20%3E%20span%2Ebn%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Efl%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Ech%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Est%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Eco%20%7B%20color%3A%20%23888888%3B%20font%2Dstyle%3A%20italic%3B%20%7D%20%0Acode%20%3E%20span%2Eot%20%7B%20color%3A%20%23007020%3B%20%7D%20%0Acode%20%3E%20span%2Eal%20%7B%20color%3A%20%23ff0000%3B%20font%2Dweight%3A%20bold%3B%20%7D%20%0Acode%20%3E%20span%2Efu%20%7B%20color%3A%20%23900%3B%20font%2Dweight%3A%20bold%3B%20%7D%20%0Acode%20%3E%20span%2Eer%20%7B%20color%3A%20%23a61717%3B%20background%2Dcolor%3A%20%23e3d2d2%3B%20%7D%20%0A" type="text/css" />




</head>

<body>




<h1 class="title toc-ignore">regulaTER: an R package for analysis of transposable elements in accessible regions</h1>



<p>. This repo is currently under review. Citation information will be provided as soon as our work is accepted.</p>
<div id="what-is-this-package-used-for" class="section level3">
<h3>What is this package used for?</h3>
<p>This package is used for analysis of repeat elements in the genome, more specifically transposable elements, and their association with accessible chromatin and gene expression regulatory elements. As input, it requires DNA accessibility data, such as that produced by ATAC-seq, and analyzed by an appropriate pipeline, a RepeatMasker file including information for genomic repeats, BED files describing the gene-related contexts of all genomic regions, and a list of differentially expressed genes in condition of choice compared to controls and their genomic coordinates.</p>
<div id="what-are-the-dependencies-for-regulater" class="section level4">
<h4>What are the dependencies for regulaTER ?</h4>
<ol style="list-style-type: decimal">
<li><a href="https://www.r-project.org/">R</a> version should be version 3.5+</li>
<li>While using R programming, we recommend the use of <a href="https://www.rstudio.com/products/rstudio/download/">Rstudio</a> for convenient use and understanding of the functions of regulaTER.</li>
<li><a href="https://bedtools.readthedocs.io/en/latest/content/installation.html">Bedtools</a> is required to be installed and in your PATH environmental variable.</li>
<li><a href="http://homer.ucsd.edu/homer/introduction/install.html">HOMER</a> is required to be installed and in your PATH environmental variable.</li>
<li><a href="https://cran.r-project.org/web/packages/devtools/readme/README.html">devtools</a> is required to install regulaTER.</li>
<li>regulaTER utilizes functions from the following R packages, which have to be installed in your R library. You may visit the following websites to install them easily:
<ul>
<li><a href="https://dplyr.tidyverse.org/">dplyr</a></li>
<li><a href="https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html">GenomicRanges</a></li>
<li><a href="https://cran.r-project.org/web/packages/biomartr/readme/README.html">biomartr</a></li>
<li><a href="https://robertamezquita.github.io/marge/index.html">marge</a></li>
</ul></li>
</ol>
</div>
</div>
<div id="how-to-install-this-r-package" class="section level3">
<h3>How to install this R package ?</h3>
<div class="sourceCode" id="cb1"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb1-1"><a href="#cb1-1" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb1-2"><a href="#cb1-2" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(devtools)</span>
<span id="cb1-3"><a href="#cb1-3" aria-hidden="true" tabindex="-1"></a>devtools<span class="sc">::</span><span class="fu">install_github</span>(<span class="st">&quot;karakulahg/regulaTER&quot;</span>)</span></code></pre></div>
</div>
<div id="how-does-it-work" class="section level3">
<h3>How does it work?</h3>
<ol style="list-style-type: decimal">
<li>Install and load the following libraries:</li>
</ol>
<div class="sourceCode" id="cb2"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a><span class="fu"></span></span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(devtools)  </span>
<span id="cb2-3"><a href="#cb2-3" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(regulaTER)</span>
<span id="cb2-4"><a href="#cb2-4" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(biomartr)</span>
<span id="cb2-5"><a href="#cb2-5" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(GenomicRanges)</span>
<span id="cb2-6"><a href="#cb2-6" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(dplyr)</span>
<span id="cb2-7"><a href="#cb2-7" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(valr)</span>
<span id="cb2-8"><a href="#cb2-8" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(biomaRt)</span>
<span id="cb2-9"><a href="#cb2-9" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(marge)</span>
<span id="cb2-10"><a href="#cb2-10" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(ChIPseeker)</span>
<span id="cb2-11"><a href="#cb2-11" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb2-12"><a href="#cb2-12" aria-hidden="true" tabindex="-1"></a><span class="fu">sessionInfo</span>()</span>
<span id="cb2-13"><a href="#cb2-13" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb2-14"><a href="#cb2-14" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb2-15"><a href="#cb2-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; R version 3.6.0 (2019-04-26)</span></span>
<span id="cb2-16"><a href="#cb2-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Platform: x86_64-redhat-linux-gnu (64-bit)</span></span>
<span id="cb2-17"><a href="#cb2-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Running under: CentOS Linux 7 (Core)</span></span>
<span id="cb2-18"><a href="#cb2-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb2-19"><a href="#cb2-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Matrix products: default</span></span>
<span id="cb2-20"><a href="#cb2-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so</span></span>
<span id="cb2-21"><a href="#cb2-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb2-22"><a href="#cb2-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; locale:</span></span>
<span id="cb2-23"><a href="#cb2-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              </span></span>
<span id="cb2-24"><a href="#cb2-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [3] LC_TIME=tr_TR.UTF-8        LC_COLLATE=en_US.UTF-8    </span></span>
<span id="cb2-25"><a href="#cb2-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [5] LC_MONETARY=tr_TR.UTF-8    LC_MESSAGES=en_US.UTF-8   </span></span>
<span id="cb2-26"><a href="#cb2-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [7] LC_PAPER=tr_TR.UTF-8       LC_NAME=C                 </span></span>
<span id="cb2-27"><a href="#cb2-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [9] LC_ADDRESS=C               LC_TELEPHONE=C            </span></span>
<span id="cb2-28"><a href="#cb2-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [11] LC_MEASUREMENT=tr_TR.UTF-8 LC_IDENTIFICATION=C       </span></span>
<span id="cb2-29"><a href="#cb2-29" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb2-30"><a href="#cb2-30" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; attached base packages:</span></span>
<span id="cb2-31"><a href="#cb2-31" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] parallel  stats4    stats     graphics  grDevices utils     datasets </span></span>
<span id="cb2-32"><a href="#cb2-32" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [8] methods   base     </span></span>
<span id="cb2-33"><a href="#cb2-33" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb2-34"><a href="#cb2-34" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; other attached packages:</span></span>
<span id="cb2-35"><a href="#cb2-35" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [1] ChIPseeker_1.22.1    marge_0.0.4.9999     biomaRt_2.42.1      </span></span>
<span id="cb2-36"><a href="#cb2-36" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [4] valr_0.6.2           dplyr_1.0.5          GenomicRanges_1.38.0</span></span>
<span id="cb2-37"><a href="#cb2-37" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [7] GenomeInfoDb_1.22.1  IRanges_2.20.2       S4Vectors_0.24.4    </span></span>
<span id="cb2-38"><a href="#cb2-38" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [10] BiocGenerics_0.32.0  biomartr_0.9.2       regulaTER_0.1.0     </span></span>
<span id="cb2-39"><a href="#cb2-39" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [13] devtools_2.3.2       usethis_2.0.1       </span></span>
<span id="cb2-40"><a href="#cb2-40" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb2-41"><a href="#cb2-41" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; loaded via a namespace (and not attached):</span></span>
<span id="cb2-42"><a href="#cb2-42" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1] backports_1.2.1                        </span></span>
<span id="cb2-43"><a href="#cb2-43" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2] fastmatch_1.1-0                        </span></span>
<span id="cb2-44"><a href="#cb2-44" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3] BiocFileCache_1.10.2                   </span></span>
<span id="cb2-45"><a href="#cb2-45" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4] plyr_1.8.6                             </span></span>
<span id="cb2-46"><a href="#cb2-46" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5] igraph_1.2.6                           </span></span>
<span id="cb2-47"><a href="#cb2-47" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6] lazyeval_0.2.2                         </span></span>
<span id="cb2-48"><a href="#cb2-48" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [7] splines_3.6.0                          </span></span>
<span id="cb2-49"><a href="#cb2-49" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [8] BiocParallel_1.20.1                    </span></span>
<span id="cb2-50"><a href="#cb2-50" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [9] ggplot2_3.3.3                          </span></span>
<span id="cb2-51"><a href="#cb2-51" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [10] urltools_1.7.3                         </span></span>
<span id="cb2-52"><a href="#cb2-52" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [11] digest_0.6.27                          </span></span>
<span id="cb2-53"><a href="#cb2-53" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [12] htmltools_0.5.1.1                      </span></span>
<span id="cb2-54"><a href="#cb2-54" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [13] GOSemSim_2.12.1                        </span></span>
<span id="cb2-55"><a href="#cb2-55" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [14] viridis_0.5.1                          </span></span>
<span id="cb2-56"><a href="#cb2-56" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [15] GO.db_3.10.0                           </span></span>
<span id="cb2-57"><a href="#cb2-57" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [16] fansi_0.4.2                            </span></span>
<span id="cb2-58"><a href="#cb2-58" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [17] magrittr_2.0.1                         </span></span>
<span id="cb2-59"><a href="#cb2-59" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [18] memoise_2.0.0                          </span></span>
<span id="cb2-60"><a href="#cb2-60" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [19] remotes_2.2.0                          </span></span>
<span id="cb2-61"><a href="#cb2-61" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [20] Biostrings_2.54.0                      </span></span>
<span id="cb2-62"><a href="#cb2-62" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [21] readr_1.4.0                            </span></span>
<span id="cb2-63"><a href="#cb2-63" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [22] graphlayouts_0.7.1                     </span></span>
<span id="cb2-64"><a href="#cb2-64" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [23] matrixStats_0.58.0                     </span></span>
<span id="cb2-65"><a href="#cb2-65" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [24] askpass_1.1                            </span></span>
<span id="cb2-66"><a href="#cb2-66" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [25] enrichplot_1.6.1                       </span></span>
<span id="cb2-67"><a href="#cb2-67" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [26] prettyunits_1.1.1                      </span></span>
<span id="cb2-68"><a href="#cb2-68" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [27] colorspace_2.0-0                       </span></span>
<span id="cb2-69"><a href="#cb2-69" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [28] blob_1.2.1                             </span></span>
<span id="cb2-70"><a href="#cb2-70" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [29] rappdirs_0.3.3                         </span></span>
<span id="cb2-71"><a href="#cb2-71" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [30] ggrepel_0.9.1                          </span></span>
<span id="cb2-72"><a href="#cb2-72" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [31] xfun_0.21                              </span></span>
<span id="cb2-73"><a href="#cb2-73" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [32] callr_3.5.1                            </span></span>
<span id="cb2-74"><a href="#cb2-74" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [33] crayon_1.4.1                           </span></span>
<span id="cb2-75"><a href="#cb2-75" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [34] RCurl_1.98-1.2                         </span></span>
<span id="cb2-76"><a href="#cb2-76" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [35] jsonlite_1.7.2                         </span></span>
<span id="cb2-77"><a href="#cb2-77" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [36] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2</span></span>
<span id="cb2-78"><a href="#cb2-78" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [37] glue_1.4.2                             </span></span>
<span id="cb2-79"><a href="#cb2-79" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [38] polyclip_1.10-0                        </span></span>
<span id="cb2-80"><a href="#cb2-80" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [39] gtable_0.3.0                           </span></span>
<span id="cb2-81"><a href="#cb2-81" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [40] zlibbioc_1.32.0                        </span></span>
<span id="cb2-82"><a href="#cb2-82" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [41] XVector_0.26.0                         </span></span>
<span id="cb2-83"><a href="#cb2-83" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [42] DelayedArray_0.12.3                    </span></span>
<span id="cb2-84"><a href="#cb2-84" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [43] pkgbuild_1.2.0                         </span></span>
<span id="cb2-85"><a href="#cb2-85" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [44] scales_1.1.1                           </span></span>
<span id="cb2-86"><a href="#cb2-86" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [45] DOSE_3.12.0                            </span></span>
<span id="cb2-87"><a href="#cb2-87" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [46] DBI_1.1.1                              </span></span>
<span id="cb2-88"><a href="#cb2-88" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [47] Rcpp_1.0.6                             </span></span>
<span id="cb2-89"><a href="#cb2-89" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [48] plotrix_3.8-1                          </span></span>
<span id="cb2-90"><a href="#cb2-90" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [49] viridisLite_0.3.0                      </span></span>
<span id="cb2-91"><a href="#cb2-91" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [50] progress_1.2.2                         </span></span>
<span id="cb2-92"><a href="#cb2-92" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [51] gridGraphics_0.5-1                     </span></span>
<span id="cb2-93"><a href="#cb2-93" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [52] bit_4.0.4                              </span></span>
<span id="cb2-94"><a href="#cb2-94" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [53] europepmc_0.4                          </span></span>
<span id="cb2-95"><a href="#cb2-95" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [54] httr_1.4.2                             </span></span>
<span id="cb2-96"><a href="#cb2-96" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [55] fgsea_1.12.0                           </span></span>
<span id="cb2-97"><a href="#cb2-97" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [56] gplots_3.1.1                           </span></span>
<span id="cb2-98"><a href="#cb2-98" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [57] RColorBrewer_1.1-2                     </span></span>
<span id="cb2-99"><a href="#cb2-99" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [58] ellipsis_0.3.1                         </span></span>
<span id="cb2-100"><a href="#cb2-100" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [59] pkgconfig_2.0.3                        </span></span>
<span id="cb2-101"><a href="#cb2-101" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [60] XML_3.99-0.3                           </span></span>
<span id="cb2-102"><a href="#cb2-102" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [61] farver_2.1.0                           </span></span>
<span id="cb2-103"><a href="#cb2-103" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [62] dbplyr_2.1.0                           </span></span>
<span id="cb2-104"><a href="#cb2-104" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [63] utf8_1.1.4                             </span></span>
<span id="cb2-105"><a href="#cb2-105" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [64] ggplotify_0.0.5                        </span></span>
<span id="cb2-106"><a href="#cb2-106" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [65] tidyselect_1.1.0                       </span></span>
<span id="cb2-107"><a href="#cb2-107" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [66] rlang_0.4.11                           </span></span>
<span id="cb2-108"><a href="#cb2-108" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [67] reshape2_1.4.4                         </span></span>
<span id="cb2-109"><a href="#cb2-109" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [68] AnnotationDbi_1.48.0                   </span></span>
<span id="cb2-110"><a href="#cb2-110" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [69] munsell_0.5.0                          </span></span>
<span id="cb2-111"><a href="#cb2-111" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [70] tools_3.6.0                            </span></span>
<span id="cb2-112"><a href="#cb2-112" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [71] cachem_1.0.4                           </span></span>
<span id="cb2-113"><a href="#cb2-113" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [72] cli_2.3.1                              </span></span>
<span id="cb2-114"><a href="#cb2-114" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [73] generics_0.1.0                         </span></span>
<span id="cb2-115"><a href="#cb2-115" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [74] RSQLite_2.2.3                          </span></span>
<span id="cb2-116"><a href="#cb2-116" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [75] ggridges_0.5.3                         </span></span>
<span id="cb2-117"><a href="#cb2-117" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [76] broom_0.7.5                            </span></span>
<span id="cb2-118"><a href="#cb2-118" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [77] evaluate_0.14                          </span></span>
<span id="cb2-119"><a href="#cb2-119" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [78] stringr_1.4.0                          </span></span>
<span id="cb2-120"><a href="#cb2-120" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [79] fastmap_1.1.0                          </span></span>
<span id="cb2-121"><a href="#cb2-121" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [80] yaml_2.2.1                             </span></span>
<span id="cb2-122"><a href="#cb2-122" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [81] processx_3.4.5                         </span></span>
<span id="cb2-123"><a href="#cb2-123" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [82] knitr_1.31                             </span></span>
<span id="cb2-124"><a href="#cb2-124" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [83] bit64_4.0.5                            </span></span>
<span id="cb2-125"><a href="#cb2-125" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [84] fs_1.5.0                               </span></span>
<span id="cb2-126"><a href="#cb2-126" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [85] tidygraph_1.2.0                        </span></span>
<span id="cb2-127"><a href="#cb2-127" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [86] caTools_1.18.1                         </span></span>
<span id="cb2-128"><a href="#cb2-128" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [87] purrr_0.3.4                            </span></span>
<span id="cb2-129"><a href="#cb2-129" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [88] ggraph_2.0.5                           </span></span>
<span id="cb2-130"><a href="#cb2-130" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [89] xml2_1.3.2                             </span></span>
<span id="cb2-131"><a href="#cb2-131" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [90] DO.db_2.9                              </span></span>
<span id="cb2-132"><a href="#cb2-132" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [91] compiler_3.6.0                         </span></span>
<span id="cb2-133"><a href="#cb2-133" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [92] curl_4.3                               </span></span>
<span id="cb2-134"><a href="#cb2-134" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [93] testthat_3.0.2                         </span></span>
<span id="cb2-135"><a href="#cb2-135" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [94] tibble_3.1.0                           </span></span>
<span id="cb2-136"><a href="#cb2-136" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [95] tweenr_1.0.1                           </span></span>
<span id="cb2-137"><a href="#cb2-137" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [96] stringi_1.5.3                          </span></span>
<span id="cb2-138"><a href="#cb2-138" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [97] ps_1.6.0                               </span></span>
<span id="cb2-139"><a href="#cb2-139" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [98] GenomicFeatures_1.38.2                 </span></span>
<span id="cb2-140"><a href="#cb2-140" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [99] desc_1.3.0                             </span></span>
<span id="cb2-141"><a href="#cb2-141" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [100] lattice_0.20-38                        </span></span>
<span id="cb2-142"><a href="#cb2-142" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [101] Matrix_1.2-17                          </span></span>
<span id="cb2-143"><a href="#cb2-143" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [102] vctrs_0.3.6                            </span></span>
<span id="cb2-144"><a href="#cb2-144" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [103] pillar_1.5.1                           </span></span>
<span id="cb2-145"><a href="#cb2-145" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [104] lifecycle_1.0.0                        </span></span>
<span id="cb2-146"><a href="#cb2-146" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [105] BiocManager_1.30.10                    </span></span>
<span id="cb2-147"><a href="#cb2-147" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [106] triebeard_0.3.0                        </span></span>
<span id="cb2-148"><a href="#cb2-148" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [107] data.table_1.14.0                      </span></span>
<span id="cb2-149"><a href="#cb2-149" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [108] cowplot_1.1.1                          </span></span>
<span id="cb2-150"><a href="#cb2-150" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [109] bitops_1.0-6                           </span></span>
<span id="cb2-151"><a href="#cb2-151" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [110] rtracklayer_1.46.0                     </span></span>
<span id="cb2-152"><a href="#cb2-152" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [111] qvalue_2.18.0                          </span></span>
<span id="cb2-153"><a href="#cb2-153" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [112] R6_2.5.0                               </span></span>
<span id="cb2-154"><a href="#cb2-154" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [113] KernSmooth_2.23-15                     </span></span>
<span id="cb2-155"><a href="#cb2-155" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [114] gridExtra_2.3                          </span></span>
<span id="cb2-156"><a href="#cb2-156" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [115] sessioninfo_1.1.1                      </span></span>
<span id="cb2-157"><a href="#cb2-157" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [116] gtools_3.8.2                           </span></span>
<span id="cb2-158"><a href="#cb2-158" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [117] boot_1.3-22                            </span></span>
<span id="cb2-159"><a href="#cb2-159" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [118] MASS_7.3-51.4                          </span></span>
<span id="cb2-160"><a href="#cb2-160" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [119] assertthat_0.2.1                       </span></span>
<span id="cb2-161"><a href="#cb2-161" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [120] pkgload_1.2.0                          </span></span>
<span id="cb2-162"><a href="#cb2-162" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [121] SummarizedExperiment_1.16.1            </span></span>
<span id="cb2-163"><a href="#cb2-163" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [122] openssl_1.4.3                          </span></span>
<span id="cb2-164"><a href="#cb2-164" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [123] rprojroot_2.0.2                        </span></span>
<span id="cb2-165"><a href="#cb2-165" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [124] withr_2.4.1                            </span></span>
<span id="cb2-166"><a href="#cb2-166" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [125] GenomicAlignments_1.22.1               </span></span>
<span id="cb2-167"><a href="#cb2-167" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [126] Rsamtools_2.2.3                        </span></span>
<span id="cb2-168"><a href="#cb2-168" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [127] GenomeInfoDbData_1.2.2                 </span></span>
<span id="cb2-169"><a href="#cb2-169" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [128] hms_1.0.0                              </span></span>
<span id="cb2-170"><a href="#cb2-170" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [129] grid_3.6.0                             </span></span>
<span id="cb2-171"><a href="#cb2-171" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [130] tidyr_1.1.3                            </span></span>
<span id="cb2-172"><a href="#cb2-172" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [131] rvcheck_0.1.8                          </span></span>
<span id="cb2-173"><a href="#cb2-173" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [132] rmarkdown_2.7                          </span></span>
<span id="cb2-174"><a href="#cb2-174" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [133] ggforce_0.3.3                          </span></span>
<span id="cb2-175"><a href="#cb2-175" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [134] Biobase_2.46.0</span></span></code></pre></div>
<ol start="2" style="list-style-type: decimal">
<li>Read the repeat track information for your genome of interest in RepeatMasker file format, as obtained from the <a href="https://www.girinst.org/repbase/">RepBase database</a>, as well as your genomic sequencing peak information, generated by MACS2 and annotated by ChIPseeker, in one of narrowPeak, broadPeak, or summits formats. For details on the format of annotated BED files, see “ChIPseeker Annotations Compatible with ShufflePeaks and TEAR Functions” below these instructions.</li>
</ol>
<div class="sourceCode" id="cb3"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a><span class="fu"></span><span class="at"></span> <span class="sc"></span><span class="dv"></span></span>
<span id="cb3-2"><a href="#cb3-2" aria-hidden="true" tabindex="-1"></a>raw.rmsk<span class="ot">&lt;-</span>biomartr<span class="sc">::</span><span class="fu">read_rm</span>(<span class="st">&quot;../hg38/hg38.fa.out.gz&quot;</span>)</span>
<span id="cb3-3"><a href="#cb3-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 80924 out of 5622516 rows ~ 0.014% were removed from the imported RepeatMasker file, because they contained &#39;NA&#39; values in either &#39;qry_start&#39; or &#39;qry_end&#39;.</span></span>
<span id="cb3-4"><a href="#cb3-4" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(raw.rmsk)</span>
<span id="cb3-5"><a href="#cb3-5" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; # A tibble: 6 x 15</span></span>
<span id="cb3-6"><a href="#cb3-6" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   sw_score perc_div perc_del perc_insert qry_id qry_start qry_end qry_left   </span></span>
<span id="cb3-7"><a href="#cb3-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   &lt;chr&gt;    &lt;chr&gt;    &lt;chr&gt;    &lt;chr&gt;       &lt;chr&gt;      &lt;int&gt;   &lt;int&gt; &lt;chr&gt;      </span></span>
<span id="cb3-8"><a href="#cb3-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 1 463      1.3      0.6      1.7         chr1       10001   10468 (248945954)</span></span>
<span id="cb3-9"><a href="#cb3-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2 4005     11.3     21.5     1.3         chr1       10469   11447 (248944975)</span></span>
<span id="cb3-10"><a href="#cb3-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3 535      21.2     15.9     3.1         chr1       11485   11676 (248944746)</span></span>
<span id="cb3-11"><a href="#cb3-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 4 263      29.4     1.9      1.0         chr1       11678   11780 (248944642)</span></span>
<span id="cb3-12"><a href="#cb3-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 5 309      23.0     3.7      0.0         chr1       15265   15355 (248941067)</span></span>
<span id="cb3-13"><a href="#cb3-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 6 18       23.2     0.0      2.0         chr1       15798   15849 (248940573)</span></span>
<span id="cb3-14"><a href="#cb3-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; # … with 7 more variables: matching_repeat &lt;chr&gt;, repeat_id &lt;chr&gt;,</span></span>
<span id="cb3-15"><a href="#cb3-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; #   matching_class &lt;chr&gt;, no_bp_in_complement &lt;chr&gt;, in_repeat_start &lt;chr&gt;,</span></span>
<span id="cb3-16"><a href="#cb3-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; #   in_repeat_end &lt;chr&gt;, qry_width &lt;int&gt;</span></span>
<span id="cb3-17"><a href="#cb3-17" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb3-18"><a href="#cb3-18" aria-hidden="true" tabindex="-1"></a>peak <span class="ot">&lt;-</span><span class="fu">readPeakFile</span>(<span class="st">&quot;../Annotated_FinalPeakList.narrowPeak.filt.gz.tsv&quot;</span>)</span>
<span id="cb3-19"><a href="#cb3-19" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb3-20"><a href="#cb3-20" aria-hidden="true" tabindex="-1"></a><span class="co"># explicitly identify the column stating the distance of the summit nucleotide from the peak start nucleotide (for narrowPeak files only)</span></span>
<span id="cb3-21"><a href="#cb3-21" aria-hidden="true" tabindex="-1"></a><span class="fu">names</span>(<span class="fu">mcols</span>(peak))[<span class="dv">8</span>] <span class="ot">&lt;-</span> <span class="st">&quot;summit&quot;</span></span>
<span id="cb3-22"><a href="#cb3-22" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(peak)</span>
<span id="cb3-23"><a href="#cb3-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; GRanges object with 6 ranges and 16 metadata columns:</span></span>
<span id="cb3-24"><a href="#cb3-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;       seqnames              ranges strand |     width   strand          V4</span></span>
<span id="cb3-25"><a href="#cb3-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;          &lt;Rle&gt;           &lt;IRanges&gt;  &lt;Rle&gt; | &lt;integer&gt; &lt;factor&gt;    &lt;factor&gt;</span></span>
<span id="cb3-26"><a href="#cb3-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1]    chr10 100009817-100010338      * |       523        * Peak_291566</span></span>
<span id="cb3-27"><a href="#cb3-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2]    chr10 100009817-100010338      * |       523        * Peak_663918</span></span>
<span id="cb3-28"><a href="#cb3-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3]    chr10 100014102-100014314      * |       214        *  Peak_37465</span></span>
<span id="cb3-29"><a href="#cb3-29" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4]    chr10 100047315-100047555      * |       242        * Peak_534187</span></span>
<span id="cb3-30"><a href="#cb3-30" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5]    chr10 100057263-100058004      * |       743        *     Peak_69</span></span>
<span id="cb3-31"><a href="#cb3-31" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6]    chr10 100095216-100095748      * |       534        *   Peak_9036</span></span>
<span id="cb3-32"><a href="#cb3-32" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;              V5        V7          V8          V9    summit</span></span>
<span id="cb3-33"><a href="#cb3-33" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;       &lt;integer&gt; &lt;numeric&gt;   &lt;numeric&gt;   &lt;numeric&gt; &lt;integer&gt;</span></span>
<span id="cb3-34"><a href="#cb3-34" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1]       177   4.00186     17.7356    15.91056       312</span></span>
<span id="cb3-35"><a href="#cb3-35" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2]        95   2.98321     9.58673     8.16916       209</span></span>
<span id="cb3-36"><a href="#cb3-36" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3]       441    6.5485    44.16472     41.4122       108</span></span>
<span id="cb3-37"><a href="#cb3-37" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4]       119   3.09789    11.90047    10.36341       124</span></span>
<span id="cb3-38"><a href="#cb3-38" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5]      1000  39.50494 35842.92578 35837.49609       284</span></span>
<span id="cb3-39"><a href="#cb3-39" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6]       707   8.65858    70.70222    67.35577       138</span></span>
<span id="cb3-40"><a href="#cb3-40" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;                                           annotation   geneChr geneStart</span></span>
<span id="cb3-41"><a href="#cb3-41" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;                                             &lt;factor&gt; &lt;integer&gt; &lt;integer&gt;</span></span>
<span id="cb3-42"><a href="#cb3-42" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1]                               Promoter (&lt;=1kb)        10  99875577</span></span>
<span id="cb3-43"><a href="#cb3-43" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2]                               Promoter (&lt;=1kb)        10  99875577</span></span>
<span id="cb3-44"><a href="#cb3-44" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3]                              Distal Intergenic        10  99875577</span></span>
<span id="cb3-45"><a href="#cb3-45" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4] Intron (ENST00000370418.8/1369, intron 8 of 8)        10 100042193</span></span>
<span id="cb3-46"><a href="#cb3-46" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5] Intron (ENST00000370418.8/1369, intron 5 of 8)        10 100042193</span></span>
<span id="cb3-47"><a href="#cb3-47" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6]                              Distal Intergenic        10 100042193</span></span>
<span id="cb3-48"><a href="#cb3-48" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;         geneEnd geneLength geneStrand    geneId distanceToTSS</span></span>
<span id="cb3-49"><a href="#cb3-49" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;       &lt;integer&gt;  &lt;integer&gt;  &lt;integer&gt; &lt;integer&gt;     &lt;integer&gt;</span></span>
<span id="cb3-50"><a href="#cb3-50" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1] 100009947     134371          2     23268             0</span></span>
<span id="cb3-51"><a href="#cb3-51" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2] 100009947     134371          2     23268             0</span></span>
<span id="cb3-52"><a href="#cb3-52" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3] 100009947     134371          2     23268         -4154</span></span>
<span id="cb3-53"><a href="#cb3-53" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4] 100081869      39677          2      1369         34314</span></span>
<span id="cb3-54"><a href="#cb3-54" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5] 100081869      39677          2      1369         23865</span></span>
<span id="cb3-55"><a href="#cb3-55" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6] 100081869      39677          2      1369        -13346</span></span>
<span id="cb3-56"><a href="#cb3-56" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   -------</span></span>
<span id="cb3-57"><a href="#cb3-57" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   seqinfo: 30 sequences from an unspecified genome; no seqlengths</span></span></code></pre></div>
<ol start="3" style="list-style-type: decimal">
<li>Using the Table Browser tool of UCSC Genome Browser, or another similar method, generate BED files for the following genomic regions, as compatible with ChIPseeker annotations. Promoter region is defined as 3000 bases upstream of the gene region, while Downstream region is defined as 3000 bases downstream. The genomeSizePath should point to a file with two columns, with the first column corresponding to chromosome names in order, and the second column corresponding to chromosome size. Details on how to use Table Browser to obtain the BED files is included below the instructions for regulaTER. For the UCSC human and mouse genome assemblies hg38, hg19, and mm10, the pre-generated region files are available for download via the GitHub repository <a href="https://github.com/karakulahg/regulaTER-regions">regulaTER-regions</a>.</li>
</ol>
<div class="sourceCode" id="cb4"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb4-1"><a href="#cb4-1" aria-hidden="true" tabindex="-1"></a>base_path <span class="ot">&lt;-</span> <span class="st">&quot;~/hg38_regions&quot;</span></span></code></pre></div>
<div class="sourceCode" id="cb5"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb5-1"><a href="#cb5-1" aria-hidden="true" tabindex="-1"></a>pathList <span class="ot">&lt;-</span> <span class="fu">list</span>(<span class="st">&quot;Promoter&quot;</span> <span class="ot">=</span> <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38_promoter_complement.bed&quot;</span>),</span>
<span id="cb5-2"><a href="#cb5-2" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;Exon&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38_exons_complement.bed&quot;</span>),</span>
<span id="cb5-3"><a href="#cb5-3" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;Intron&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38_introns_complement.bed&quot;</span>),</span>
<span id="cb5-4"><a href="#cb5-4" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;5UTR&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38_5prime_complement.bed&quot;</span>),</span>
<span id="cb5-5"><a href="#cb5-5" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;3UTR&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38_3prime_complement.bed&quot;</span>),</span>
<span id="cb5-6"><a href="#cb5-6" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;Downstream&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38_downstream_complement.bed&quot;</span>),</span>
<span id="cb5-7"><a href="#cb5-7" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;genomeSizePath&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/hg38.chrom.sizes&quot;</span>))</span></code></pre></div>
<ol start="4" style="list-style-type: decimal">
<li>With the previously generated peak file and repeat masker objects, as well as the list of file paths to BED files defining genomic regions and the file with chromosome sizes, calculate the fold enrichment and p-values of repeat elements associated with input peak regions. If a broadPeak file is used, the format should be changed as “broad”, and minoverlap should be given an integer value, defining how many minimum bases should a peak and a repeat element overlap before they are considered associated. Higher shuffle numbers give better results, but increase operation time.</li>
</ol>
<div class="sourceCode" id="cb6"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb6-1"><a href="#cb6-1" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-2"><a href="#cb6-2" aria-hidden="true" tabindex="-1"></a>enrich.peak <span class="ot">&lt;-</span> <span class="fu">TEAR</span>(<span class="at">inputPeakFile =</span> peak, <span class="at">pathList =</span> pathList, <span class="at">numberOfShuffle =</span> <span class="dv">2</span>, <span class="at">repeatMaskerFile =</span> raw.rmsk, <span class="at">format=</span><span class="st">&quot;narrow&quot;</span>, <span class="at">minoverlap=</span>0L, <span class="at">alternative =</span> <span class="st">&quot;greater&quot;</span>, <span class="at">minobserved =</span> <span class="dv">10</span>)</span>
<span id="cb6-3"><a href="#cb6-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] &quot;shuffle  1  - sh time : 8.12676572799683&quot;</span></span>
<span id="cb6-4"><a href="#cb6-4" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] &quot;shuffle  2  - sh time : 8.22793889045715&quot;</span></span>
<span id="cb6-5"><a href="#cb6-5" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-6"><a href="#cb6-6" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(enrich.peak<span class="sc">$</span>RepeatName)</span>
<span id="cb6-7"><a href="#cb6-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     RepeatName   rmsk observed expected TrueMean      p.value p.adjust.value</span></span>
<span id="cb6-8"><a href="#cb6-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 324      L1P4b    153       23        0     0.25 0.000000e+00   0.000000e+00</span></span>
<span id="cb6-9"><a href="#cb6-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 938      MLT1K  17999      236       88    88.25 3.280682e-39   9.886381e-38</span></span>
<span id="cb6-10"><a href="#cb6-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 362        L2b  93700      817      511   511.00 5.477844e-36   1.607315e-34</span></span>
<span id="cb6-11"><a href="#cb6-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 335  L1PA15-16   1195       77       14    14.00 2.300437e-32   6.576890e-31</span></span>
<span id="cb6-12"><a href="#cb6-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 520     LTR41B    868       53        6     6.25 3.163887e-32   8.819336e-31</span></span>
<span id="cb6-13"><a href="#cb6-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 361        L2a 170012     1529     1119  1119.00 1.327943e-31   3.611357e-30</span></span>
<span id="cb6-14"><a href="#cb6-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     obsOnTrueMean</span></span>
<span id="cb6-15"><a href="#cb6-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 324     92.000000</span></span>
<span id="cb6-16"><a href="#cb6-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 938      2.674221</span></span>
<span id="cb6-17"><a href="#cb6-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 362      1.598826</span></span>
<span id="cb6-18"><a href="#cb6-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 335      5.500000</span></span>
<span id="cb6-19"><a href="#cb6-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 520      8.480000</span></span>
<span id="cb6-20"><a href="#cb6-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 361      1.366399</span></span></code></pre></div>
<ol start="5" style="list-style-type: decimal">
<li>Using the TEAR results, the peak file and repeat masker objects used as the input for TEAR, and a user provided data frame of differentially expressed genes and their genomic intervals (such as those obtained from ensembl’s BioMart or UCSC’s Genome Browser), the structure of which is described under the function documentation, calculate the fold enrichment and p-values of accessible region enriched repeat elements located in input promoter regions. Higher shuffle numbers give better results, but increase operation time. Distance value defines desired length of promoter region from transcription start site. Longer distances cover more distal regulatory elements, but reduce precision of results.</li>
</ol>
<div class="sourceCode" id="cb7"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb7-1"><a href="#cb7-1" aria-hidden="true" tabindex="-1"></a>input <span class="ot">&lt;-</span> <span class="st">&quot;../DE_genenames.txt&quot;</span></span>
<span id="cb7-2"><a href="#cb7-2" aria-hidden="true" tabindex="-1"></a>dataset <span class="ot">&lt;-</span> <span class="st">&quot;hsapiens_gene_ensembl&quot;</span></span>
<span id="cb7-3"><a href="#cb7-3" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-4"><a href="#cb7-4" aria-hidden="true" tabindex="-1"></a>genes <span class="ot">&lt;-</span> regulaTER<span class="sc">::</span><span class="fu">getInterval</span>(input, dataset)</span>
<span id="cb7-5"><a href="#cb7-5" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(genes)</span>
<span id="cb7-6"><a href="#cb7-6" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;            geneID geneName seqnames     start       end strand</span></span>
<span id="cb7-7"><a href="#cb7-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 1 ENSG00000000003   TSPAN6        X 100627108 100639991      -</span></span>
<span id="cb7-8"><a href="#cb7-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2 ENSG00000000419     DPM1       20  50934867  50959140      -</span></span>
<span id="cb7-9"><a href="#cb7-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3 ENSG00000000457    SCYL3     chr1 169849631 169894267      -</span></span>
<span id="cb7-10"><a href="#cb7-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 4 ENSG00000000460 C1orf112     chr1 169662007 169854080      +</span></span>
<span id="cb7-11"><a href="#cb7-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 5 ENSG00000001084     GCLC     chr6  53497341  53616970      -</span></span>
<span id="cb7-12"><a href="#cb7-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 6 ENSG00000001167     NFYA     chr6  41072974  41102403      +</span></span>
<span id="cb7-13"><a href="#cb7-13" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-14"><a href="#cb7-14" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-15"><a href="#cb7-15" aria-hidden="true" tabindex="-1"></a>IdDEGRepeats <span class="ot">&lt;-</span> <span class="fu">DATE</span>(<span class="at">enrichTEARResult =</span> enrich.peak, <span class="at">peaks =</span> peak, <span class="at">rmsk =</span> raw.rmsk, <span class="at">genes =</span> genes, <span class="at">alternative =</span> <span class="st">&quot;greater&quot;</span> ,<span class="at">numberOfShuffle =</span> <span class="dv">2</span>,<span class="at">minobserved =</span> <span class="dv">10</span>, <span class="at">distance =</span> <span class="dv">100000</span>)</span>
<span id="cb7-16"><a href="#cb7-16" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(IdDEGRepeats)</span>
<span id="cb7-17"><a href="#cb7-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;    RepeatName rmsk observed expected p.value p.adjust.value</span></span>
<span id="cb7-18"><a href="#cb7-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 9      AluSc5 4259        1        0      NA             NA</span></span>
<span id="cb7-19"><a href="#cb7-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 30      AluYc 6889        1        0      NA             NA</span></span>
<span id="cb7-20"><a href="#cb7-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 46     AluYk2 6900        1        0      NA             NA</span></span>
<span id="cb7-21"><a href="#cb7-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 51   AmnSINE2  671        1        0      NA             NA</span></span>
<span id="cb7-22"><a href="#cb7-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 53   Arthur1A 1183        1        0      NA             NA</span></span>
<span id="cb7-23"><a href="#cb7-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 55   Arthur1C  727        1        0      NA             NA</span></span></code></pre></div>
<ol start="6" style="list-style-type: decimal">
<li>Using the output of DATE, as well as the peak, repeat masker, and genes object used in the same function, identify which genomic motifs are enriched in promoter associated TEAR regions. This function requires a local installation of HOMER, available on the <a href="http://homer.ucsd.edu/homer/">HOMER website</a>.</li>
</ol>
<div class="sourceCode" id="cb8"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb8-1"><a href="#cb8-1" aria-hidden="true" tabindex="-1"></a><span class="fu">FindMotifs</span>(<span class="at">df =</span> result100, <span class="at">repeatMaskerFile =</span> raw.rmsk, <span class="at">peak =</span> peak, <span class="at">distance =</span> <span class="dv">100000</span>, <span class="at">genes =</span> genes, <span class="at">genome =</span> <span class="st">&quot;hg38&quot;</span>, <span class="at">outDir =</span> <span class="st">&quot;../test/&quot;</span>, <span class="at">homerPath =</span> <span class="st">&quot;~/Downloads/Tools/homer/&quot;</span>, <span class="at">type =</span> <span class="st">&quot;linkedRepeats&quot;</span>, <span class="at">numberOfMotifs =</span> T)</span></code></pre></div>
</div>
<div id="using-the-ucsc-table-browser-to-obtain-genomic-regions-compatible-with-chipseeker-regions" class="section level3">
<h3>Using the UCSC Table Browser to Obtain Genomic Regions Compatible with ChIPseeker Regions</h3>
<p>ChIPseeker regions are divided into seven categories. In order of priority, these are promoter, coding exon, intron, 5’ UTR, 3’ UTR, downstream regions, and intergenic regions. Overlaps are categorized based on the order of priority when annotating genomic regions. As the shuffling of peak regions are performed in the same category of region as they are originally found in, the intervals for these regions must be supplied by the user. If you are working with a genome with known gene annotations available on the UCSC Genome Browser, you can obtain these files using the Table Browser tool.</p>
<p>Under genome and assembly, select the organism and genome assembly of choice, also to be used with ChIPseeker to annotate the peak regions. For region, select “genome”, and for output format, choose “BED - browser extensible data”. Fill the output file field to download the resulting file to your system.</p>
<p>On the next screen, you can further select the regions your file will have. For promoter and downstream regions, select “Upstream by” or “Downstream by”, respectively, in both cases using 3000 bases as the region length. Intergenic regions need not be downloaded, as any region not falling under another category are included in intergenic shuffles.</p>
<p>In the current version of regulaTER, the complement intervals of the regions must be provided for each object in the list. They use the same format as the obtained BED files, but cover all genomic regions not found in the category. The user can generate these files using another genome arithmetic tool, such as bedtools complement (-i <BED> -g <genome>), available on the <a href="https://bedtools.readthedocs.io/en/latest/index.html">bedtools Suite</a>.</p>
<p>The user will also require the full lengths of the chromosomes included in the annotation. These are available for each genome found in UCSC Genome Browser under Downloads / Genome Data, with the name “[assembly].chrom.sizes”, under the link named “Genome sequence files and select annotations (2bit, GTF, GC-content, etc)”.</p>
<p>If the user is using annotations not available on the UCSC Genome Browser, they can manually generate these files using the tool of their choice, as long as gene and exon annotations are already available in the BED format.</p>
</div>
<div id="chipseeker-annotations-compatible-with-shufflepeaks-and-tear-functions" class="section level3">
<h3>ChIPseeker Annotations Compatible with ShufflePeaks and TEAR Functions</h3>
<p>As of MACS2 (version X) and ChIPseeker (version X), a narrowPeak format BED file generated by MACS2 and annotated by ChIPseeker has 19 columns. The first eleven columns of the file correspond to the original BED file, with the remaining columns added by ChIPseeker. Only the first give columns, the summit column and the annotation fields are used during regulaTER’s shufflePeak function. seqnames identifies the chromosome the peak is located in, start and end columns identify the 5’ and 3’ of the interval on the sense strand, width identifies the length of the interval, strand identifies the strand information of the peak, if applicable. The “summit” column (by default the 13th column of the BED file) identifies the distance of the summit nucleotide from the interval start site. The “annotation” column identifies the peak’s relationship to nearby genes, and can be one of Promoter, Exon, Intron, 5’ UTR, 3’ UTR, Downstream, and Distal Intergenic, with further annotation specifying distance from TSS or TES, as well as which exon or intron, if applicable. If the user desires to supply their own annotated interval regions, the minimum required columns for the GRanges object used in regulaTER functions are named “seqnames”, “ranges” (“start” and “end”), “strand”, “summit” and “annotation”.</p>
<p>BED files in broadPeak or summits have one fewer column than those in the narrowPeak format, as they lack the summit column.</p>
</div>



<!-- code folding -->


<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

</body>
</html>
