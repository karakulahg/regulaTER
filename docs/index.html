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
<p>For the analysis of narrowPeak peak files generated by MACS2, follow the steps outlined below.</p>
<ol style="list-style-type: decimal">
<li>Install and load the following libraries:</li>
</ol>
<div class="sourceCode" id="cb2"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb2-1"><a href="#cb2-1" aria-hidden="true" tabindex="-1"></a><span class="fu"></span></span>
<span id="cb2-2"><a href="#cb2-2" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb2-3"><a href="#cb2-3" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(regulaTER)</span>
<span id="cb2-4"><a href="#cb2-4" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(biomartr)</span>
<span id="cb2-5"><a href="#cb2-5" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(GenomicRanges)</span>
<span id="cb2-6"><a href="#cb2-6" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(dplyr)</span>
<span id="cb2-7"><a href="#cb2-7" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(biomaRt)</span>
<span id="cb2-8"><a href="#cb2-8" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(marge)</span>
<span id="cb2-9"><a href="#cb2-9" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(ChIPseeker)</span>
<span id="cb2-10"><a href="#cb2-10" aria-hidden="true" tabindex="-1"></a>  </span>
<span id="cb2-11"><a href="#cb2-11" aria-hidden="true" tabindex="-1"></a></span></code></pre></div>
<ol start="2" style="list-style-type: decimal">
<li>Read the repeat track information for your genome of interest in RepeatMasker file format, as obtained from the <a href="https://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html">RepeatMasker database</a>, as well as your genomic sequencing peak information, generated by MACS2 and annotated by ChIPseeker, in one of narrowPeak, broadPeak, or summits formats. For details on the format of annotated BED files, see “ChIPseeker Annotations Compatible with ShufflePeaks and TEAR Functions” below these instructions.</li>
</ol>
<div class="sourceCode" id="cb3"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a><span class="fu"></span><span class="at"></span> <span class="sc"></span><span class="dv"></span></span>
<span id="cb3-2"><a href="#cb3-2" aria-hidden="true" tabindex="-1"></a>raw.rmsk<span class="ot">&lt;-</span>biomartr<span class="sc">::</span><span class="fu">read_rm</span>(<span class="st">&quot;~/hg38/hg38.fa.out.gz&quot;</span>)</span>
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
<span id="cb3-18"><a href="#cb3-18" aria-hidden="true" tabindex="-1"></a>peak <span class="ot">&lt;-</span><span class="fu">readPeakFile</span>(<span class="st">&quot;Annotated_FinalPeakList.narrowPeak.filt.tsv&quot;</span>)</span>
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
<li>With the previously generated peak file and repeat masker objects, as well as the list of file paths to BED files defining genomic regions and the file with chromosome sizes, calculate the fold enrichment and p-values of repeat elements associated with input peak regions (Transposable Elements in Accessible Regions - TEAR). If a broadPeak file is used, the format should be changed as “broad”, and minoverlap should be given an integer value, defining how many minimum bases should a peak and a repeat element overlap before they are considered associated. Higher shuffle numbers give better results, but increase operation time.</li>
</ol>
<div class="sourceCode" id="cb6"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb6-1"><a href="#cb6-1" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-2"><a href="#cb6-2" aria-hidden="true" tabindex="-1"></a>enrich.peak <span class="ot">&lt;-</span> <span class="fu">TEAR</span>(<span class="at">inputPeakFile =</span> peak, <span class="at">pathList =</span> pathList, <span class="at">numberOfShuffle =</span> <span class="dv">2</span>, <span class="at">repeatMaskerFile =</span> raw.rmsk, <span class="at">format=</span><span class="st">&quot;narrow&quot;</span>, <span class="at">minoverlap=</span>0L, <span class="at">alternative =</span> <span class="st">&quot;greater&quot;</span>, <span class="at">minobserved =</span> <span class="dv">10</span>)</span>
<span id="cb6-3"><a href="#cb6-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] &quot;shuffle  1  - sh time : 8.35915327072144&quot;</span></span>
<span id="cb6-4"><a href="#cb6-4" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] &quot;shuffle  2  - sh time : 8.34707999229431&quot;</span></span>
<span id="cb6-5"><a href="#cb6-5" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-6"><a href="#cb6-6" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(enrich.peak<span class="sc">$</span>RepeatName)<span class="fu"></span><span class="sc"></span><span class="sc"></span></span>
<span id="cb6-7"><a href="#cb6-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     RepeatName  rmsk observed expected TrueMean      p.value p.adjust.value</span></span>
<span id="cb6-8"><a href="#cb6-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 233     L1M2a1   116       11        0     0.25 0.000000e+00   0.000000e+00</span></span>
<span id="cb6-9"><a href="#cb6-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 520     LTR41B   868       53        5     5.25 5.076533e-36   1.825914e-34</span></span>
<span id="cb6-10"><a href="#cb6-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 938      MLT1K 17999      236       96    96.25 1.026175e-33   3.575578e-32</span></span>
<span id="cb6-11"><a href="#cb6-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 362        L2b 93700      817      523   523.25 6.265120e-33   2.116851e-31</span></span>
<span id="cb6-12"><a href="#cb6-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 491      LTR33  9406      147       47    47.00 1.210416e-31   3.969452e-30</span></span>
<span id="cb6-13"><a href="#cb6-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 418    LTR16A2  1849       54        7     7.25 1.079059e-29   3.437573e-28</span></span>
<span id="cb6-14"><a href="#cb6-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     obsOnTrueMean</span></span>
<span id="cb6-15"><a href="#cb6-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 233     44.000000</span></span>
<span id="cb6-16"><a href="#cb6-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 520     10.095238</span></span>
<span id="cb6-17"><a href="#cb6-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 938      2.451948</span></span>
<span id="cb6-18"><a href="#cb6-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 362      1.561395</span></span>
<span id="cb6-19"><a href="#cb6-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 491      3.127660</span></span>
<span id="cb6-20"><a href="#cb6-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 418      7.448276</span></span></code></pre></div>
<ol start="5" style="list-style-type: decimal">
<li>Using the results of the TEAR function, the peak file and repeat masker objects used as the input for TEAR, and a user provided data frame of differentially expressed genes and their genomic intervals (such as those obtained from ensembl’s BioMart or UCSC’s Genome Browser), the structure of which is described under the function documentation, calculate the fold enrichment and p-values of accessible region enriched repeat elements located in input promoter regions (DEG Associated Transposable Elements - DATE). Higher shuffle numbers give better results, but increase operation time. Distance value defines desired length of promoter region from transcription start site. Longer distances cover more distal regulatory elements, but reduce precision of results.</li>
</ol>
<div class="sourceCode" id="cb7"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb7-1"><a href="#cb7-1" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-2"><a href="#cb7-2" aria-hidden="true" tabindex="-1"></a><span class="co"># use getInterval to generate interval ranges from gene list</span></span>
<span id="cb7-3"><a href="#cb7-3" aria-hidden="true" tabindex="-1"></a>input <span class="ot">&lt;-</span> <span class="st">&quot;DE_genenames.txt&quot;</span></span>
<span id="cb7-4"><a href="#cb7-4" aria-hidden="true" tabindex="-1"></a>dataset <span class="ot">&lt;-</span> <span class="st">&quot;hsapiens_gene_ensembl&quot;</span></span>
<span id="cb7-5"><a href="#cb7-5" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-6"><a href="#cb7-6" aria-hidden="true" tabindex="-1"></a>genes <span class="ot">&lt;-</span> regulaTER<span class="sc">::</span><span class="fu">getInterval</span>(input, dataset)</span>
<span id="cb7-7"><a href="#cb7-7" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(genes)</span>
<span id="cb7-8"><a href="#cb7-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;            geneID geneName seqnames     start       end strand</span></span>
<span id="cb7-9"><a href="#cb7-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 1 ENSG00000000003   TSPAN6        X 100627108 100639991      -</span></span>
<span id="cb7-10"><a href="#cb7-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2 ENSG00000000419     DPM1       20  50934867  50959140      -</span></span>
<span id="cb7-11"><a href="#cb7-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3 ENSG00000000457    SCYL3     chr1 169849631 169894267      -</span></span>
<span id="cb7-12"><a href="#cb7-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 4 ENSG00000000460 C1orf112     chr1 169662007 169854080      +</span></span>
<span id="cb7-13"><a href="#cb7-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 5 ENSG00000001084     GCLC     chr6  53497341  53616970      -</span></span>
<span id="cb7-14"><a href="#cb7-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 6 ENSG00000001167     NFYA     chr6  41072974  41102403      +</span></span>
<span id="cb7-15"><a href="#cb7-15" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-16"><a href="#cb7-16" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-17"><a href="#cb7-17" aria-hidden="true" tabindex="-1"></a>IdDEGRepeats <span class="ot">&lt;-</span> <span class="fu">DATE</span>(<span class="at">enrichTEARResult =</span> enrich.peak, <span class="at">peaks =</span> peak, <span class="at">rmsk =</span> raw.rmsk, <span class="at">genes =</span> genes, <span class="at">alternative =</span> <span class="st">&quot;greater&quot;</span> ,<span class="at">numberOfShuffle =</span> <span class="dv">2</span>, <span class="at">minobserved =</span> <span class="dv">10</span>, <span class="at">distance =</span> <span class="dv">100000</span>)</span>
<span id="cb7-18"><a href="#cb7-18" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb7-19"><a href="#cb7-19" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(IdDEGRepeats)<span class="fu"></span><span class="sc"></span></span>
<span id="cb7-20"><a href="#cb7-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;           RepeatName   rmsk observed expected      p.value p.adjust.value</span></span>
<span id="cb7-21"><a href="#cb7-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 355               L2  55157      125       82 6.007228e-06   7.442288e-05</span></span>
<span id="cb7-22"><a href="#cb7-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 123  ERV3-16A3_I-int   6654       35       15 7.084875e-06   8.680918e-05</span></span>
<span id="cb7-23"><a href="#cb7-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 938            MLT1K  17999       57       32 4.163017e-05   5.045396e-04</span></span>
<span id="cb7-24"><a href="#cb7-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 361              L2a 170012      331      267 8.543413e-05   1.024291e-03</span></span>
<span id="cb7-25"><a href="#cb7-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 884             MIRb 228856      264      209 1.392829e-04   1.652133e-03</span></span>
<span id="cb7-26"><a href="#cb7-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 1000           SVA_D   1232       18        7 3.450966e-04   4.050344e-03</span></span></code></pre></div>
<ol start="6" style="list-style-type: decimal">
<li>Using the output of the DATE function, as well as the peak, repeat masker, and genes object used in the same function, identify which genomic motifs are enriched in promoter associated TEAR regions. This function requires a local installation of HOMER, available on the <a href="http://homer.ucsd.edu/homer/">HOMER website</a>.</li>
</ol>
<div class="sourceCode" id="cb8"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb8-1"><a href="#cb8-1" aria-hidden="true" tabindex="-1"></a><span class="fu">FindMotifs</span>(<span class="at">df =</span> IdDEGRepeats, <span class="at">repeatMaskerFile =</span> raw.rmsk, <span class="at">peak =</span> peak, <span class="at">distance =</span> <span class="dv">100000</span>, <span class="at">genes =</span> genes, <span class="at">genome =</span> <span class="st">&quot;hg38&quot;</span>, <span class="at">outDir =</span> <span class="st">&quot;../test/&quot;</span>, <span class="at">homerPath =</span> <span class="st">&quot;~/homer/&quot;</span>, <span class="at">type =</span> <span class="st">&quot;linkedRepeats&quot;</span>, <span class="at">topRepeats =</span> T)</span></code></pre></div>
<p>Alternately, the repeatName object in the output of TEAR can be provided instead of the output of DATE. In this case, the distance and genes parameters can be omitted.</p>
<div class="sourceCode" id="cb9"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb9-1"><a href="#cb9-1" aria-hidden="true" tabindex="-1"></a><span class="fu">FindMotifs</span>(<span class="at">df =</span> enrich.peak<span class="sc">$</span>RepeatName, <span class="at">repeatMaskerFile =</span> raw.rmsk, <span class="at">peak =</span> peak, <span class="at">genome =</span> <span class="st">&quot;hg38&quot;</span>, <span class="at">outDir =</span> <span class="st">&quot;../test/&quot;</span>, <span class="at">homerPath =</span> <span class="st">&quot;~/Downloads/Tools/homer/&quot;</span>, <span class="at">type =</span> <span class="st">&quot;enrichPeak&quot;</span>, <span class="at">topRepeats =</span> T)</span></code></pre></div>
</div>
<div id="analysis-of-broadpeak-and-summits-files" class="section level3">
<h3>Analysis of broadPeak and summits files</h3>
<p>For broadPeak files, many of the above steps can be used as is. As broad peak analysis does not identify summit nucleotide, broadPeak files lack the summit column. Therefore, the summit column is not used for the analysis, and does not have to be identified. In addition, the minimum amount of overlapping bases between the peaks and the repeat regions has to be provided as a long integer. It is recommended that the value is no larger than half the length of the shortest repeat of interest, or half the read length of the sequencing library used to identify peaks, whichever is smaller.</p>
<div class="sourceCode" id="cb10"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb10-1"><a href="#cb10-1" aria-hidden="true" tabindex="-1"></a><span class="fu"></span><span class="at"></span> <span class="sc"></span><span class="dv"></span></span>
<span id="cb10-2"><a href="#cb10-2" aria-hidden="true" tabindex="-1"></a>raw.rmsk<span class="ot">&lt;-</span>biomartr<span class="sc">::</span><span class="fu">read_rm</span>(<span class="st">&quot;~/mm10/mm10.fa.out&quot;</span>)</span>
<span id="cb10-3"><a href="#cb10-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 76798 out of 5294587 rows ~ 0.015% were removed from the imported RepeatMasker file, because they contained &#39;NA&#39; values in either &#39;qry_start&#39; or &#39;qry_end&#39;.</span></span>
<span id="cb10-4"><a href="#cb10-4" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(raw.rmsk)</span>
<span id="cb10-5"><a href="#cb10-5" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; # A tibble: 6 x 15</span></span>
<span id="cb10-6"><a href="#cb10-6" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   sw_score perc_div perc_del perc_insert qry_id qry_start qry_end qry_left   </span></span>
<span id="cb10-7"><a href="#cb10-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   &lt;chr&gt;    &lt;chr&gt;    &lt;chr&gt;    &lt;chr&gt;       &lt;chr&gt;      &lt;int&gt;   &lt;int&gt; &lt;chr&gt;      </span></span>
<span id="cb10-8"><a href="#cb10-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 1 27       0.0      0.0      0.0         chr1     3000098 3000123 (192471848)</span></span>
<span id="cb10-9"><a href="#cb10-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 2 2853     18.6     1.9      6.8         chr1     3003148 3004054 (192467917)</span></span>
<span id="cb10-10"><a href="#cb10-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 3 3936     19.9     2.1      1.4         chr1     3004041 3004206 (192467765)</span></span>
<span id="cb10-11"><a href="#cb10-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 4 43       6.5      3.1      0.0         chr1     3004207 3004270 (192467701)</span></span>
<span id="cb10-12"><a href="#cb10-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 5 3936     19.9     2.1      1.4         chr1     3004271 3005001 (192466970)</span></span>
<span id="cb10-13"><a href="#cb10-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 6 1392     22.3     4.4      6.1         chr1     3005002 3005441 (192466530)</span></span>
<span id="cb10-14"><a href="#cb10-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; # … with 7 more variables: matching_repeat &lt;chr&gt;, repeat_id &lt;chr&gt;,</span></span>
<span id="cb10-15"><a href="#cb10-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; #   matching_class &lt;chr&gt;, no_bp_in_complement &lt;chr&gt;, in_repeat_start &lt;chr&gt;,</span></span>
<span id="cb10-16"><a href="#cb10-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; #   in_repeat_end &lt;chr&gt;, qry_width &lt;int&gt;</span></span>
<span id="cb10-17"><a href="#cb10-17" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb10-18"><a href="#cb10-18" aria-hidden="true" tabindex="-1"></a>peak <span class="ot">&lt;-</span><span class="fu">readPeakFile</span>(<span class="st">&quot;Annotated_Vehicle_K27Ac_IP-S1_R1_peaks.broadPeak.tsv&quot;</span>)</span>
<span id="cb10-19"><a href="#cb10-19" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb10-20"><a href="#cb10-20" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(peak)</span>
<span id="cb10-21"><a href="#cb10-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; GRanges object with 6 ranges and 15 metadata columns:</span></span>
<span id="cb10-22"><a href="#cb10-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;       seqnames          ranges strand |     width   strand</span></span>
<span id="cb10-23"><a href="#cb10-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;          &lt;Rle&gt;       &lt;IRanges&gt;  &lt;Rle&gt; | &lt;integer&gt; &lt;factor&gt;</span></span>
<span id="cb10-24"><a href="#cb10-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1]     chr1 3004745-3004891      * |       148        *</span></span>
<span id="cb10-25"><a href="#cb10-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2]     chr1 3008558-3008695      * |       139        *</span></span>
<span id="cb10-26"><a href="#cb10-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3]     chr1 3010320-3010443      * |       125        *</span></span>
<span id="cb10-27"><a href="#cb10-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4]     chr1 3014326-3014523      * |       199        *</span></span>
<span id="cb10-28"><a href="#cb10-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5]     chr1 3027534-3027576      * |        44        *</span></span>
<span id="cb10-29"><a href="#cb10-29" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6]     chr1 3028186-3028609      * |       425        *</span></span>
<span id="cb10-30"><a href="#cb10-30" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;                                  V4        V5        V7        V8        V9</span></span>
<span id="cb10-31"><a href="#cb10-31" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;                            &lt;factor&gt; &lt;integer&gt; &lt;numeric&gt; &lt;numeric&gt; &lt;numeric&gt;</span></span>
<span id="cb10-32"><a href="#cb10-32" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1] Vehicle_K27Ac_IP-S1_R1_peak_1        18   1.74877   1.87941         0</span></span>
<span id="cb10-33"><a href="#cb10-33" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2] Vehicle_K27Ac_IP-S1_R1_peak_2        21   1.81713   2.12406         0</span></span>
<span id="cb10-34"><a href="#cb10-34" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3] Vehicle_K27Ac_IP-S1_R1_peak_3        21   1.81713   2.12406         0</span></span>
<span id="cb10-35"><a href="#cb10-35" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4] Vehicle_K27Ac_IP-S1_R1_peak_4        13   1.54366   1.36369         0</span></span>
<span id="cb10-36"><a href="#cb10-36" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5] Vehicle_K27Ac_IP-S1_R1_peak_5        12   1.67066   1.27923         0</span></span>
<span id="cb10-37"><a href="#cb10-37" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6] Vehicle_K27Ac_IP-S1_R1_peak_6        13   1.53074   1.33756         0</span></span>
<span id="cb10-38"><a href="#cb10-38" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;              annotation   geneChr geneStart   geneEnd geneLength geneStrand</span></span>
<span id="cb10-39"><a href="#cb10-39" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;                &lt;factor&gt; &lt;integer&gt; &lt;integer&gt; &lt;integer&gt;  &lt;integer&gt;  &lt;integer&gt;</span></span>
<span id="cb10-40"><a href="#cb10-40" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1] Distal Intergenic         1   3214482   3671498     457017          2</span></span>
<span id="cb10-41"><a href="#cb10-41" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2] Distal Intergenic         1   3214482   3671498     457017          2</span></span>
<span id="cb10-42"><a href="#cb10-42" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3] Distal Intergenic         1   3214482   3671498     457017          2</span></span>
<span id="cb10-43"><a href="#cb10-43" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4] Distal Intergenic         1   3214482   3671498     457017          2</span></span>
<span id="cb10-44"><a href="#cb10-44" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5] Distal Intergenic         1   3214482   3671498     457017          2</span></span>
<span id="cb10-45"><a href="#cb10-45" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6] Distal Intergenic         1   3214482   3671498     457017          2</span></span>
<span id="cb10-46"><a href="#cb10-46" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;          geneId distanceToTSS</span></span>
<span id="cb10-47"><a href="#cb10-47" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;       &lt;integer&gt;     &lt;numeric&gt;</span></span>
<span id="cb10-48"><a href="#cb10-48" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1]    497097        666607</span></span>
<span id="cb10-49"><a href="#cb10-49" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2]    497097        662803</span></span>
<span id="cb10-50"><a href="#cb10-50" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3]    497097        661055</span></span>
<span id="cb10-51"><a href="#cb10-51" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4]    497097        656975</span></span>
<span id="cb10-52"><a href="#cb10-52" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5]    497097        643922</span></span>
<span id="cb10-53"><a href="#cb10-53" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6]    497097        642889</span></span>
<span id="cb10-54"><a href="#cb10-54" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   -------</span></span>
<span id="cb10-55"><a href="#cb10-55" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   seqinfo: 21 sequences from an unspecified genome; no seqlengths</span></span>
<span id="cb10-56"><a href="#cb10-56" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb10-57"><a href="#cb10-57" aria-hidden="true" tabindex="-1"></a>base_path <span class="ot">&lt;-</span> <span class="st">&quot;~/mm10_regions&quot;</span></span>
<span id="cb10-58"><a href="#cb10-58" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb10-59"><a href="#cb10-59" aria-hidden="true" tabindex="-1"></a>pathList <span class="ot">&lt;-</span> <span class="fu">list</span>(<span class="st">&quot;Promoter&quot;</span> <span class="ot">=</span> <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_promoter_complement.bed&quot;</span>),</span>
<span id="cb10-60"><a href="#cb10-60" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;Exon&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_exons_complement.bed&quot;</span>),</span>
<span id="cb10-61"><a href="#cb10-61" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;Intron&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_introns_complement.bed&quot;</span>),</span>
<span id="cb10-62"><a href="#cb10-62" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;5UTR&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_5prime_complement.bed&quot;</span>),</span>
<span id="cb10-63"><a href="#cb10-63" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;3UTR&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_3prime_complement.bed&quot;</span>),</span>
<span id="cb10-64"><a href="#cb10-64" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;Downstream&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_downstream_complement.bed&quot;</span>),</span>
<span id="cb10-65"><a href="#cb10-65" aria-hidden="true" tabindex="-1"></a>                 <span class="st">&quot;genomeSizePath&quot;</span> <span class="ot">=</span>  <span class="fu">paste0</span>(base_path, <span class="st">&quot;/mm10_chrom.sizes&quot;</span>))</span>
<span id="cb10-66"><a href="#cb10-66" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb10-67"><a href="#cb10-67" aria-hidden="true" tabindex="-1"></a>enrich.peak <span class="ot">&lt;-</span> <span class="fu">TEAR</span>(<span class="at">inputPeakFile =</span> peak, <span class="at">pathList =</span> pathList, <span class="at">numberOfShuffle =</span> <span class="dv">2</span>, <span class="at">repeatMaskerFile =</span> raw.rmsk, <span class="at">format=</span><span class="st">&quot;broad&quot;</span>, <span class="at">minoverlap=</span>10L, <span class="at">alternative =</span> <span class="st">&quot;greater&quot;</span>, <span class="at">minobserved =</span> <span class="dv">10</span>)</span>
<span id="cb10-68"><a href="#cb10-68" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] &quot;shuffle  1  - sh time : 51.6609251499176&quot;</span></span>
<span id="cb10-69"><a href="#cb10-69" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] &quot;shuffle  2  - sh time : 51.8320786952972&quot;</span></span>
<span id="cb10-70"><a href="#cb10-70" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb10-71"><a href="#cb10-71" aria-hidden="true" tabindex="-1"></a><span class="fu">head</span>(enrich.peak<span class="sc">$</span>RepeatName)<span class="fu"></span><span class="sc"></span><span class="sc"></span></span>
<span id="cb10-72"><a href="#cb10-72" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     RepeatName  rmsk observed expected TrueMean p.value p.adjust.value</span></span>
<span id="cb10-73"><a href="#cb10-73" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 195    L1_Mur2 17465     3379     1634  1634.25       0              0</span></span>
<span id="cb10-74"><a href="#cb10-74" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 196    L1_Mur3 25665     5562     2861  2861.25       0              0</span></span>
<span id="cb10-75"><a href="#cb10-75" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 451        Lx7 39622     7365     4118  4117.75       0              0</span></span>
<span id="cb10-76"><a href="#cb10-76" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 454        Lx9 49209     8434     5135  5134.75       0              0</span></span>
<span id="cb10-77"><a href="#cb10-77" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 855 RLTR17B_Mm 29320     5376     2813  2813.25       0              0</span></span>
<span id="cb10-78"><a href="#cb10-78" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 881   RLTR20A4 15354     3313     1270  1270.00       0              0</span></span>
<span id="cb10-79"><a href="#cb10-79" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;     obsOnTrueMean</span></span>
<span id="cb10-80"><a href="#cb10-80" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 195      2.067615</span></span>
<span id="cb10-81"><a href="#cb10-81" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 196      1.943906</span></span>
<span id="cb10-82"><a href="#cb10-82" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 451      1.788598</span></span>
<span id="cb10-83"><a href="#cb10-83" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 454      1.642534</span></span>
<span id="cb10-84"><a href="#cb10-84" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 855      1.910957</span></span>
<span id="cb10-85"><a href="#cb10-85" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; 881      2.608661</span></span></code></pre></div>
<p>Files in the summits format only define the summit nucleotide for each identified peak, and use the same analysis steps as for narrowPeak files, with the exception that the summit column does not have to be identified.</p>
</div>
<div id="using-the-ucsc-table-browser-to-obtain-genomic-regions-compatible-with-chipseeker-regions" class="section level3">
<h3>Using the UCSC Table Browser to Obtain Genomic Regions Compatible with ChIPseeker Regions</h3>
<p>ChIPseeker regions are divided into seven categories. In order of priority, these are promoter, coding exon, intron, 5’ UTR, 3’ UTR, downstream regions, and intergenic regions. Overlaps are categorized based on the order of priority when annotating genomic regions. As the shuffling of peak regions are performed in the same category of region as they are originally found in, the intervals for these regions must be supplied by the user. If you are working with a genome with known gene annotations available on the UCSC Genome Browser, you can obtain these files using the Table Browser tool.</p>
<p>Under genome and assembly, select the organism and genome assembly of choice, also to be used with ChIPseeker to annotate the peak regions. For region, select “genome”, and for output format, choose “BED - browser extensible data”. Fill the output file field to download the resulting file to your system.</p>
<p>On the next screen, you can further select the regions your file will have. For promoter and downstream regions, select “Upstream by” or “Downstream by”, respectively, in both cases using 3000 bases as the region length. Intergenic regions need not be downloaded, as any region not falling under another category are included in intergenic shuffles.</p>
<p>In the current version of regulaTER, the complement intervals of the regions must be provided for each object in the list. They use the same format as the obtained BED files, but cover all genomic regions not found in the category. The user can generate these files using another genome arithmetic tool, such as bedtools complement (-i BED -g genome), available on the <a href="https://bedtools.readthedocs.io/en/latest/index.html">bedtools Suite</a>.</p>
<p>The user will also require the full lengths of the chromosomes included in the annotation. These are available for each genome found in UCSC Genome Browser under Downloads / Genome Data, with the name “[assembly].chrom.sizes”, under the link named “Genome sequence files and select annotations (2bit, GTF, GC-content, etc)”.</p>
<p>If the user is using annotations not available on the UCSC Genome Browser, they can manually generate these files using the tool of their choice, as long as gene and exon annotations are already available in the BED format.</p>
</div>
<div id="chipseeker-annotations-compatible-with-shufflepeaks-and-tear-functions" class="section level3">
<h3>ChIPseeker Annotations Compatible with ShufflePeaks and TEAR Functions</h3>
<p>As of MACS2 (version 2.1.2) and ChIPseeker (version 1.22.1), a narrowPeak format BED file generated by MACS2 and annotated by ChIPseeker has 19 columns. The first eleven columns of the file correspond to the original BED file, with the remaining columns added by ChIPseeker. Only the first give columns, the summit column and the annotation fields are used during regulaTER’s shufflePeak function. seqnames identifies the chromosome the peak is located in, start and end columns identify the 5’ and 3’ of the interval on the sense strand, width identifies the length of the interval, strand identifies the strand information of the peak, if applicable. The “summit” column (by default the 13th column of the BED file) identifies the distance of the summit nucleotide from the interval start site. The “annotation” column identifies the peak’s relationship to nearby genes, and can be one of Promoter, Exon, Intron, 5’ UTR, 3’ UTR, Downstream, and Distal Intergenic, with further annotation specifying distance from TSS or TES, as well as which exon or intron, if applicable. If the user desires to supply their own annotated interval regions, the minimum required columns for the GRanges object used in regulaTER functions are named “seqnames”, “ranges” (“start” and “end”), “strand”, “summit” and “annotation”.</p>
<p>BED files in broadPeak or summits have one fewer column than those in the narrowPeak format, as they lack the summit column.</p>
</div>
<div id="session-information" class="section level3">
<h3>Session Information</h3>
<div class="sourceCode" id="cb11"><pre class="sourceCode r"><code class="sourceCode r"><span id="cb11-1"><a href="#cb11-1" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb11-2"><a href="#cb11-2" aria-hidden="true" tabindex="-1"></a><span class="fu">sessionInfo</span>()</span>
<span id="cb11-3"><a href="#cb11-3" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; R version 3.6.0 (2019-04-26)</span></span>
<span id="cb11-4"><a href="#cb11-4" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Platform: x86_64-redhat-linux-gnu (64-bit)</span></span>
<span id="cb11-5"><a href="#cb11-5" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Running under: CentOS Linux 7 (Core)</span></span>
<span id="cb11-6"><a href="#cb11-6" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb11-7"><a href="#cb11-7" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; Matrix products: default</span></span>
<span id="cb11-8"><a href="#cb11-8" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so</span></span>
<span id="cb11-9"><a href="#cb11-9" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb11-10"><a href="#cb11-10" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; locale:</span></span>
<span id="cb11-11"><a href="#cb11-11" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              </span></span>
<span id="cb11-12"><a href="#cb11-12" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [3] LC_TIME=tr_TR.UTF-8        LC_COLLATE=en_US.UTF-8    </span></span>
<span id="cb11-13"><a href="#cb11-13" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [5] LC_MONETARY=tr_TR.UTF-8    LC_MESSAGES=en_US.UTF-8   </span></span>
<span id="cb11-14"><a href="#cb11-14" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [7] LC_PAPER=tr_TR.UTF-8       LC_NAME=C                 </span></span>
<span id="cb11-15"><a href="#cb11-15" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [9] LC_ADDRESS=C               LC_TELEPHONE=C            </span></span>
<span id="cb11-16"><a href="#cb11-16" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [11] LC_MEASUREMENT=tr_TR.UTF-8 LC_IDENTIFICATION=C       </span></span>
<span id="cb11-17"><a href="#cb11-17" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb11-18"><a href="#cb11-18" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; attached base packages:</span></span>
<span id="cb11-19"><a href="#cb11-19" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [1] parallel  stats4    stats     graphics  grDevices utils     datasets </span></span>
<span id="cb11-20"><a href="#cb11-20" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [8] methods   base     </span></span>
<span id="cb11-21"><a href="#cb11-21" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb11-22"><a href="#cb11-22" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; other attached packages:</span></span>
<span id="cb11-23"><a href="#cb11-23" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [1] ChIPseeker_1.22.1    marge_0.0.4.9999     biomaRt_2.42.1      </span></span>
<span id="cb11-24"><a href="#cb11-24" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [4] dplyr_1.0.5          GenomicRanges_1.38.0 GenomeInfoDb_1.22.1 </span></span>
<span id="cb11-25"><a href="#cb11-25" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [7] IRanges_2.20.2       S4Vectors_0.24.4     BiocGenerics_0.32.0 </span></span>
<span id="cb11-26"><a href="#cb11-26" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [10] biomartr_0.9.2       regulaTER_0.1.0     </span></span>
<span id="cb11-27"><a href="#cb11-27" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; </span></span>
<span id="cb11-28"><a href="#cb11-28" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; loaded via a namespace (and not attached):</span></span>
<span id="cb11-29"><a href="#cb11-29" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [1] fgsea_1.12.0                           </span></span>
<span id="cb11-30"><a href="#cb11-30" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [2] colorspace_2.0-0                       </span></span>
<span id="cb11-31"><a href="#cb11-31" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [3] ggridges_0.5.3                         </span></span>
<span id="cb11-32"><a href="#cb11-32" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [4] ellipsis_0.3.1                         </span></span>
<span id="cb11-33"><a href="#cb11-33" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [5] qvalue_2.18.0                          </span></span>
<span id="cb11-34"><a href="#cb11-34" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [6] XVector_0.26.0                         </span></span>
<span id="cb11-35"><a href="#cb11-35" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [7] rstudioapi_0.13                        </span></span>
<span id="cb11-36"><a href="#cb11-36" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [8] farver_2.1.0                           </span></span>
<span id="cb11-37"><a href="#cb11-37" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;   [9] urltools_1.7.3                         </span></span>
<span id="cb11-38"><a href="#cb11-38" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [10] graphlayouts_0.7.1                     </span></span>
<span id="cb11-39"><a href="#cb11-39" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [11] ggrepel_0.9.1                          </span></span>
<span id="cb11-40"><a href="#cb11-40" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [12] bit64_4.0.5                            </span></span>
<span id="cb11-41"><a href="#cb11-41" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [13] AnnotationDbi_1.48.0                   </span></span>
<span id="cb11-42"><a href="#cb11-42" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [14] fansi_0.4.2                            </span></span>
<span id="cb11-43"><a href="#cb11-43" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [15] xml2_1.3.2                             </span></span>
<span id="cb11-44"><a href="#cb11-44" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [16] splines_3.6.0                          </span></span>
<span id="cb11-45"><a href="#cb11-45" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [17] cachem_1.0.4                           </span></span>
<span id="cb11-46"><a href="#cb11-46" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [18] GOSemSim_2.12.1                        </span></span>
<span id="cb11-47"><a href="#cb11-47" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [19] knitr_1.31                             </span></span>
<span id="cb11-48"><a href="#cb11-48" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [20] polyclip_1.10-0                        </span></span>
<span id="cb11-49"><a href="#cb11-49" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [21] jsonlite_1.7.2                         </span></span>
<span id="cb11-50"><a href="#cb11-50" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [22] Rsamtools_2.2.3                        </span></span>
<span id="cb11-51"><a href="#cb11-51" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [23] GO.db_3.10.0                           </span></span>
<span id="cb11-52"><a href="#cb11-52" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [24] dbplyr_2.1.0                           </span></span>
<span id="cb11-53"><a href="#cb11-53" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [25] ggforce_0.3.3                          </span></span>
<span id="cb11-54"><a href="#cb11-54" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [26] BiocManager_1.30.10                    </span></span>
<span id="cb11-55"><a href="#cb11-55" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [27] readr_1.4.0                            </span></span>
<span id="cb11-56"><a href="#cb11-56" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [28] compiler_3.6.0                         </span></span>
<span id="cb11-57"><a href="#cb11-57" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [29] httr_1.4.2                             </span></span>
<span id="cb11-58"><a href="#cb11-58" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [30] rvcheck_0.1.8                          </span></span>
<span id="cb11-59"><a href="#cb11-59" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [31] assertthat_0.2.1                       </span></span>
<span id="cb11-60"><a href="#cb11-60" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [32] Matrix_1.2-17                          </span></span>
<span id="cb11-61"><a href="#cb11-61" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [33] fastmap_1.1.0                          </span></span>
<span id="cb11-62"><a href="#cb11-62" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [34] lazyeval_0.2.2                         </span></span>
<span id="cb11-63"><a href="#cb11-63" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [35] cli_2.3.1                              </span></span>
<span id="cb11-64"><a href="#cb11-64" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [36] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2</span></span>
<span id="cb11-65"><a href="#cb11-65" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [37] tweenr_1.0.1                           </span></span>
<span id="cb11-66"><a href="#cb11-66" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [38] htmltools_0.5.1.1                      </span></span>
<span id="cb11-67"><a href="#cb11-67" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [39] prettyunits_1.1.1                      </span></span>
<span id="cb11-68"><a href="#cb11-68" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [40] tools_3.6.0                            </span></span>
<span id="cb11-69"><a href="#cb11-69" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [41] igraph_1.2.6                           </span></span>
<span id="cb11-70"><a href="#cb11-70" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [42] gtable_0.3.0                           </span></span>
<span id="cb11-71"><a href="#cb11-71" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [43] glue_1.4.2                             </span></span>
<span id="cb11-72"><a href="#cb11-72" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [44] GenomeInfoDbData_1.2.2                 </span></span>
<span id="cb11-73"><a href="#cb11-73" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [45] reshape2_1.4.4                         </span></span>
<span id="cb11-74"><a href="#cb11-74" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [46] DO.db_2.9                              </span></span>
<span id="cb11-75"><a href="#cb11-75" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [47] rappdirs_0.3.3                         </span></span>
<span id="cb11-76"><a href="#cb11-76" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [48] fastmatch_1.1-0                        </span></span>
<span id="cb11-77"><a href="#cb11-77" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [49] Rcpp_1.0.6                             </span></span>
<span id="cb11-78"><a href="#cb11-78" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [50] enrichplot_1.6.1                       </span></span>
<span id="cb11-79"><a href="#cb11-79" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [51] Biobase_2.46.0                         </span></span>
<span id="cb11-80"><a href="#cb11-80" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [52] vctrs_0.3.6                            </span></span>
<span id="cb11-81"><a href="#cb11-81" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [53] Biostrings_2.54.0                      </span></span>
<span id="cb11-82"><a href="#cb11-82" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [54] rtracklayer_1.46.0                     </span></span>
<span id="cb11-83"><a href="#cb11-83" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [55] ggraph_2.0.5                           </span></span>
<span id="cb11-84"><a href="#cb11-84" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [56] xfun_0.21                              </span></span>
<span id="cb11-85"><a href="#cb11-85" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [57] stringr_1.4.0                          </span></span>
<span id="cb11-86"><a href="#cb11-86" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [58] lifecycle_1.0.0                        </span></span>
<span id="cb11-87"><a href="#cb11-87" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [59] gtools_3.8.2                           </span></span>
<span id="cb11-88"><a href="#cb11-88" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [60] XML_3.99-0.3                           </span></span>
<span id="cb11-89"><a href="#cb11-89" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [61] DOSE_3.12.0                            </span></span>
<span id="cb11-90"><a href="#cb11-90" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [62] europepmc_0.4                          </span></span>
<span id="cb11-91"><a href="#cb11-91" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [63] zlibbioc_1.32.0                        </span></span>
<span id="cb11-92"><a href="#cb11-92" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [64] MASS_7.3-51.4                          </span></span>
<span id="cb11-93"><a href="#cb11-93" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [65] scales_1.1.1                           </span></span>
<span id="cb11-94"><a href="#cb11-94" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [66] tidygraph_1.2.0                        </span></span>
<span id="cb11-95"><a href="#cb11-95" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [67] hms_1.0.0                              </span></span>
<span id="cb11-96"><a href="#cb11-96" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [68] SummarizedExperiment_1.16.1            </span></span>
<span id="cb11-97"><a href="#cb11-97" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [69] RColorBrewer_1.1-2                     </span></span>
<span id="cb11-98"><a href="#cb11-98" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [70] yaml_2.2.1                             </span></span>
<span id="cb11-99"><a href="#cb11-99" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [71] curl_4.3                               </span></span>
<span id="cb11-100"><a href="#cb11-100" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [72] memoise_2.0.0                          </span></span>
<span id="cb11-101"><a href="#cb11-101" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [73] gridExtra_2.3                          </span></span>
<span id="cb11-102"><a href="#cb11-102" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [74] ggplot2_3.3.3                          </span></span>
<span id="cb11-103"><a href="#cb11-103" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [75] triebeard_0.3.0                        </span></span>
<span id="cb11-104"><a href="#cb11-104" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [76] stringi_1.5.3                          </span></span>
<span id="cb11-105"><a href="#cb11-105" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [77] RSQLite_2.2.3                          </span></span>
<span id="cb11-106"><a href="#cb11-106" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [78] plotrix_3.8-1                          </span></span>
<span id="cb11-107"><a href="#cb11-107" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [79] caTools_1.18.1                         </span></span>
<span id="cb11-108"><a href="#cb11-108" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [80] GenomicFeatures_1.38.2                 </span></span>
<span id="cb11-109"><a href="#cb11-109" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [81] boot_1.3-22                            </span></span>
<span id="cb11-110"><a href="#cb11-110" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [82] BiocParallel_1.20.1                    </span></span>
<span id="cb11-111"><a href="#cb11-111" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [83] rlang_0.4.11                           </span></span>
<span id="cb11-112"><a href="#cb11-112" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [84] pkgconfig_2.0.3                        </span></span>
<span id="cb11-113"><a href="#cb11-113" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [85] matrixStats_0.58.0                     </span></span>
<span id="cb11-114"><a href="#cb11-114" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [86] bitops_1.0-6                           </span></span>
<span id="cb11-115"><a href="#cb11-115" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [87] evaluate_0.14                          </span></span>
<span id="cb11-116"><a href="#cb11-116" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [88] lattice_0.20-38                        </span></span>
<span id="cb11-117"><a href="#cb11-117" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [89] purrr_0.3.4                            </span></span>
<span id="cb11-118"><a href="#cb11-118" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [90] GenomicAlignments_1.22.1               </span></span>
<span id="cb11-119"><a href="#cb11-119" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [91] cowplot_1.1.1                          </span></span>
<span id="cb11-120"><a href="#cb11-120" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [92] bit_4.0.4                              </span></span>
<span id="cb11-121"><a href="#cb11-121" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [93] tidyselect_1.1.0                       </span></span>
<span id="cb11-122"><a href="#cb11-122" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [94] plyr_1.8.6                             </span></span>
<span id="cb11-123"><a href="#cb11-123" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [95] magrittr_2.0.1                         </span></span>
<span id="cb11-124"><a href="#cb11-124" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [96] R6_2.5.0                               </span></span>
<span id="cb11-125"><a href="#cb11-125" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [97] gplots_3.1.1                           </span></span>
<span id="cb11-126"><a href="#cb11-126" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [98] generics_0.1.0                         </span></span>
<span id="cb11-127"><a href="#cb11-127" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt;  [99] DelayedArray_0.12.3                    </span></span>
<span id="cb11-128"><a href="#cb11-128" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [100] DBI_1.1.1                              </span></span>
<span id="cb11-129"><a href="#cb11-129" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [101] pillar_1.5.1                           </span></span>
<span id="cb11-130"><a href="#cb11-130" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [102] RCurl_1.98-1.2                         </span></span>
<span id="cb11-131"><a href="#cb11-131" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [103] tibble_3.1.0                           </span></span>
<span id="cb11-132"><a href="#cb11-132" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [104] crayon_1.4.1                           </span></span>
<span id="cb11-133"><a href="#cb11-133" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [105] KernSmooth_2.23-15                     </span></span>
<span id="cb11-134"><a href="#cb11-134" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [106] utf8_1.1.4                             </span></span>
<span id="cb11-135"><a href="#cb11-135" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [107] BiocFileCache_1.10.2                   </span></span>
<span id="cb11-136"><a href="#cb11-136" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [108] rmarkdown_2.7                          </span></span>
<span id="cb11-137"><a href="#cb11-137" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [109] viridis_0.5.1                          </span></span>
<span id="cb11-138"><a href="#cb11-138" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [110] progress_1.2.2                         </span></span>
<span id="cb11-139"><a href="#cb11-139" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [111] grid_3.6.0                             </span></span>
<span id="cb11-140"><a href="#cb11-140" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [112] data.table_1.14.0                      </span></span>
<span id="cb11-141"><a href="#cb11-141" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [113] blob_1.2.1                             </span></span>
<span id="cb11-142"><a href="#cb11-142" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [114] digest_0.6.27                          </span></span>
<span id="cb11-143"><a href="#cb11-143" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [115] tidyr_1.1.3                            </span></span>
<span id="cb11-144"><a href="#cb11-144" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [116] gridGraphics_0.5-1                     </span></span>
<span id="cb11-145"><a href="#cb11-145" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [117] openssl_1.4.3                          </span></span>
<span id="cb11-146"><a href="#cb11-146" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [118] munsell_0.5.0                          </span></span>
<span id="cb11-147"><a href="#cb11-147" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [119] ggplotify_0.0.5                        </span></span>
<span id="cb11-148"><a href="#cb11-148" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [120] viridisLite_0.3.0                      </span></span>
<span id="cb11-149"><a href="#cb11-149" aria-hidden="true" tabindex="-1"></a><span class="co">#&gt; [121] askpass_1.1</span></span></code></pre></div>
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
