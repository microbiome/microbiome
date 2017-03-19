<!DOCTYPE html>

<html xmlns="http://www.w3.org/1999/xhtml">

<head>

<meta charset="utf-8">
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="pandoc" />

<meta name="viewport" content="width=device-width, initial-scale=1">

<meta name="author" content="Sudarshan A. Shetty" />

<meta name="date" content="2017-03-05" />

<title>Format Phyloseq</title>



<style type="text/css">code{white-space: pre;}</style>
<style type="text/css">
div.sourceCode { overflow-x: auto; }
table.sourceCode, tr.sourceCode, td.lineNumbers, td.sourceCode {
  margin: 0; padding: 0; vertical-align: baseline; border: none; }
table.sourceCode { width: 100%; line-height: 100%; }
td.lineNumbers { text-align: right; padding-right: 4px; padding-left: 4px; color: #aaaaaa; border-right: 1px solid #aaaaaa; }
td.sourceCode { padding-left: 5px; }
code > span.kw { color: #007020; font-weight: bold; } /* Keyword */
code > span.dt { color: #902000; } /* DataType */
code > span.dv { color: #40a070; } /* DecVal */
code > span.bn { color: #40a070; } /* BaseN */
code > span.fl { color: #40a070; } /* Float */
code > span.ch { color: #4070a0; } /* Char */
code > span.st { color: #4070a0; } /* String */
code > span.co { color: #60a0b0; font-style: italic; } /* Comment */
code > span.ot { color: #007020; } /* Other */
code > span.al { color: #ff0000; font-weight: bold; } /* Alert */
code > span.fu { color: #06287e; } /* Function */
code > span.er { color: #ff0000; font-weight: bold; } /* Error */
code > span.wa { color: #60a0b0; font-weight: bold; font-style: italic; } /* Warning */
code > span.cn { color: #880000; } /* Constant */
code > span.sc { color: #4070a0; } /* SpecialChar */
code > span.vs { color: #4070a0; } /* VerbatimString */
code > span.ss { color: #bb6688; } /* SpecialString */
code > span.im { } /* Import */
code > span.va { color: #19177c; } /* Variable */
code > span.cf { color: #007020; font-weight: bold; } /* ControlFlow */
code > span.op { color: #666666; } /* Operator */
code > span.bu { } /* BuiltIn */
code > span.ex { } /* Extension */
code > span.pp { color: #bc7a00; } /* Preprocessor */
code > span.at { color: #7d9029; } /* Attribute */
code > span.do { color: #ba2121; font-style: italic; } /* Documentation */
code > span.an { color: #60a0b0; font-weight: bold; font-style: italic; } /* Annotation */
code > span.cv { color: #60a0b0; font-weight: bold; font-style: italic; } /* CommentVar */
code > span.in { color: #60a0b0; font-weight: bold; font-style: italic; } /* Information */
</style>



<link href="data:text/css;charset=utf-8,body%20%7B%0Abackground%2Dcolor%3A%20%23fff%3B%0Amargin%3A%201em%20auto%3B%0Amax%2Dwidth%3A%20700px%3B%0Aoverflow%3A%20visible%3B%0Apadding%2Dleft%3A%202em%3B%0Apadding%2Dright%3A%202em%3B%0Afont%2Dfamily%3A%20%22Open%20Sans%22%2C%20%22Helvetica%20Neue%22%2C%20Helvetica%2C%20Arial%2C%20sans%2Dserif%3B%0Afont%2Dsize%3A%2014px%3B%0Aline%2Dheight%3A%201%2E35%3B%0A%7D%0A%23header%20%7B%0Atext%2Dalign%3A%20center%3B%0A%7D%0A%23TOC%20%7B%0Aclear%3A%20both%3B%0Amargin%3A%200%200%2010px%2010px%3B%0Apadding%3A%204px%3B%0Awidth%3A%20400px%3B%0Aborder%3A%201px%20solid%20%23CCCCCC%3B%0Aborder%2Dradius%3A%205px%3B%0Abackground%2Dcolor%3A%20%23f6f6f6%3B%0Afont%2Dsize%3A%2013px%3B%0Aline%2Dheight%3A%201%2E3%3B%0A%7D%0A%23TOC%20%2Etoctitle%20%7B%0Afont%2Dweight%3A%20bold%3B%0Afont%2Dsize%3A%2015px%3B%0Amargin%2Dleft%3A%205px%3B%0A%7D%0A%23TOC%20ul%20%7B%0Apadding%2Dleft%3A%2040px%3B%0Amargin%2Dleft%3A%20%2D1%2E5em%3B%0Amargin%2Dtop%3A%205px%3B%0Amargin%2Dbottom%3A%205px%3B%0A%7D%0A%23TOC%20ul%20ul%20%7B%0Amargin%2Dleft%3A%20%2D2em%3B%0A%7D%0A%23TOC%20li%20%7B%0Aline%2Dheight%3A%2016px%3B%0A%7D%0Atable%20%7B%0Amargin%3A%201em%20auto%3B%0Aborder%2Dwidth%3A%201px%3B%0Aborder%2Dcolor%3A%20%23DDDDDD%3B%0Aborder%2Dstyle%3A%20outset%3B%0Aborder%2Dcollapse%3A%20collapse%3B%0A%7D%0Atable%20th%20%7B%0Aborder%2Dwidth%3A%202px%3B%0Apadding%3A%205px%3B%0Aborder%2Dstyle%3A%20inset%3B%0A%7D%0Atable%20td%20%7B%0Aborder%2Dwidth%3A%201px%3B%0Aborder%2Dstyle%3A%20inset%3B%0Aline%2Dheight%3A%2018px%3B%0Apadding%3A%205px%205px%3B%0A%7D%0Atable%2C%20table%20th%2C%20table%20td%20%7B%0Aborder%2Dleft%2Dstyle%3A%20none%3B%0Aborder%2Dright%2Dstyle%3A%20none%3B%0A%7D%0Atable%20thead%2C%20table%20tr%2Eeven%20%7B%0Abackground%2Dcolor%3A%20%23f7f7f7%3B%0A%7D%0Ap%20%7B%0Amargin%3A%200%2E5em%200%3B%0A%7D%0Ablockquote%20%7B%0Abackground%2Dcolor%3A%20%23f6f6f6%3B%0Apadding%3A%200%2E25em%200%2E75em%3B%0A%7D%0Ahr%20%7B%0Aborder%2Dstyle%3A%20solid%3B%0Aborder%3A%20none%3B%0Aborder%2Dtop%3A%201px%20solid%20%23777%3B%0Amargin%3A%2028px%200%3B%0A%7D%0Adl%20%7B%0Amargin%2Dleft%3A%200%3B%0A%7D%0Adl%20dd%20%7B%0Amargin%2Dbottom%3A%2013px%3B%0Amargin%2Dleft%3A%2013px%3B%0A%7D%0Adl%20dt%20%7B%0Afont%2Dweight%3A%20bold%3B%0A%7D%0Aul%20%7B%0Amargin%2Dtop%3A%200%3B%0A%7D%0Aul%20li%20%7B%0Alist%2Dstyle%3A%20circle%20outside%3B%0A%7D%0Aul%20ul%20%7B%0Amargin%2Dbottom%3A%200%3B%0A%7D%0Apre%2C%20code%20%7B%0Abackground%2Dcolor%3A%20%23f7f7f7%3B%0Aborder%2Dradius%3A%203px%3B%0Acolor%3A%20%23333%3B%0Awhite%2Dspace%3A%20pre%2Dwrap%3B%20%0A%7D%0Apre%20%7B%0Aborder%2Dradius%3A%203px%3B%0Amargin%3A%205px%200px%2010px%200px%3B%0Apadding%3A%2010px%3B%0A%7D%0Apre%3Anot%28%5Bclass%5D%29%20%7B%0Abackground%2Dcolor%3A%20%23f7f7f7%3B%0A%7D%0Acode%20%7B%0Afont%2Dfamily%3A%20Consolas%2C%20Monaco%2C%20%27Courier%20New%27%2C%20monospace%3B%0Afont%2Dsize%3A%2085%25%3B%0A%7D%0Ap%20%3E%20code%2C%20li%20%3E%20code%20%7B%0Apadding%3A%202px%200px%3B%0A%7D%0Adiv%2Efigure%20%7B%0Atext%2Dalign%3A%20center%3B%0A%7D%0Aimg%20%7B%0Abackground%2Dcolor%3A%20%23FFFFFF%3B%0Apadding%3A%202px%3B%0Aborder%3A%201px%20solid%20%23DDDDDD%3B%0Aborder%2Dradius%3A%203px%3B%0Aborder%3A%201px%20solid%20%23CCCCCC%3B%0Amargin%3A%200%205px%3B%0A%7D%0Ah1%20%7B%0Amargin%2Dtop%3A%200%3B%0Afont%2Dsize%3A%2035px%3B%0Aline%2Dheight%3A%2040px%3B%0A%7D%0Ah2%20%7B%0Aborder%2Dbottom%3A%204px%20solid%20%23f7f7f7%3B%0Apadding%2Dtop%3A%2010px%3B%0Apadding%2Dbottom%3A%202px%3B%0Afont%2Dsize%3A%20145%25%3B%0A%7D%0Ah3%20%7B%0Aborder%2Dbottom%3A%202px%20solid%20%23f7f7f7%3B%0Apadding%2Dtop%3A%2010px%3B%0Afont%2Dsize%3A%20120%25%3B%0A%7D%0Ah4%20%7B%0Aborder%2Dbottom%3A%201px%20solid%20%23f7f7f7%3B%0Amargin%2Dleft%3A%208px%3B%0Afont%2Dsize%3A%20105%25%3B%0A%7D%0Ah5%2C%20h6%20%7B%0Aborder%2Dbottom%3A%201px%20solid%20%23ccc%3B%0Afont%2Dsize%3A%20105%25%3B%0A%7D%0Aa%20%7B%0Acolor%3A%20%230033dd%3B%0Atext%2Ddecoration%3A%20none%3B%0A%7D%0Aa%3Ahover%20%7B%0Acolor%3A%20%236666ff%3B%20%7D%0Aa%3Avisited%20%7B%0Acolor%3A%20%23800080%3B%20%7D%0Aa%3Avisited%3Ahover%20%7B%0Acolor%3A%20%23BB00BB%3B%20%7D%0Aa%5Bhref%5E%3D%22http%3A%22%5D%20%7B%0Atext%2Ddecoration%3A%20underline%3B%20%7D%0Aa%5Bhref%5E%3D%22https%3A%22%5D%20%7B%0Atext%2Ddecoration%3A%20underline%3B%20%7D%0A%0Acode%20%3E%20span%2Ekw%20%7B%20color%3A%20%23555%3B%20font%2Dweight%3A%20bold%3B%20%7D%20%0Acode%20%3E%20span%2Edt%20%7B%20color%3A%20%23902000%3B%20%7D%20%0Acode%20%3E%20span%2Edv%20%7B%20color%3A%20%2340a070%3B%20%7D%20%0Acode%20%3E%20span%2Ebn%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Efl%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Ech%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Est%20%7B%20color%3A%20%23d14%3B%20%7D%20%0Acode%20%3E%20span%2Eco%20%7B%20color%3A%20%23888888%3B%20font%2Dstyle%3A%20italic%3B%20%7D%20%0Acode%20%3E%20span%2Eot%20%7B%20color%3A%20%23007020%3B%20%7D%20%0Acode%20%3E%20span%2Eal%20%7B%20color%3A%20%23ff0000%3B%20font%2Dweight%3A%20bold%3B%20%7D%20%0Acode%20%3E%20span%2Efu%20%7B%20color%3A%20%23900%3B%20font%2Dweight%3A%20bold%3B%20%7D%20%20code%20%3E%20span%2Eer%20%7B%20color%3A%20%23a61717%3B%20background%2Dcolor%3A%20%23e3d2d2%3B%20%7D%20%0A" rel="stylesheet" type="text/css" />

</head>

<body>




<h1 class="title toc-ignore">Format Phyloseq</h1>
<h4 class="author"><em>Sudarshan A. Shetty</em></h4>
<h4 class="date"><em>2017-03-05</em></h4>



<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - Format Phyloseq}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
<div id="formatting-the-phyloseq-object" class="section level2">
<h2>Formatting the Phyloseq Object</h2>
<p>Load <a href="Data.md">example data</a>:<br />
For this example we will use data from <a href="http://www.nature.com/articles/nmicrobiol20174">Halfvarson J., et al. Nature Microbiology, 2017</a>. It was downloaded from <a href="https://qiita.ucsd.edu/study/description/1629">Qitta</a>.</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">library</span>(microbiome)
<span class="kw">data</span>(DynamicsIBD)
p0 &lt;-<span class="st"> </span>DynamicsIBD</code></pre></div>
<div id="check-the-taxonomy" class="section level3">
<h3>Check the taxonomy</h3>
<p>We will check the taxonomy information stored in the phyloseq object.</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r"><span class="kw">library</span>(knitr)
<span class="kw">kable</span>(<span class="kw">head</span>(<span class="kw">tax_table</span>(p0)))</code></pre></div>
<table>
<thead>
<tr class="header">
<th></th>
<th align="left">Rank1</th>
<th align="left">Rank2</th>
<th align="left">Rank3</th>
<th align="left">Rank4</th>
<th align="left">Rank5</th>
<th align="left">Rank6</th>
<th align="left">Rank7</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>577110</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Clostridia</td>
<td align="left">o__Clostridiales</td>
<td align="left">f__</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td>181342</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Clostridia</td>
<td align="left">o__Clostridiales</td>
<td align="left">f__Clostridiaceae</td>
<td align="left">g__02d06</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td>581609</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Clostridia</td>
<td align="left">o__Clostridiales</td>
<td align="left">f__Ruminococcaceae</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td>4341234</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Clostridia</td>
<td align="left">o__Clostridiales</td>
<td align="left">f__Peptococcaceae</td>
<td align="left">g__Desulfotomaculum</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td>181348</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Clostridia</td>
<td align="left">o__Clostridiales</td>
<td align="left">f__Lachnospiraceae</td>
<td align="left">g__Coprococcus</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td>4467992</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Bacilli</td>
<td align="left">o__Lactobacillales</td>
<td align="left">f__Streptococcaceae</td>
<td align="left">g__Streptococcus</td>
<td align="left">s__</td>
</tr>
</tbody>
</table>
<p>It can be observed that the not all the OTUs are classified until the lowest taxonomic level (here, species level). This is especially the case with high throughput sequencing data sets. In doing OTU level testing for differential abundance, you may need information regading the specific otu number or taxonomy of the otu. This can help in easily tracing back the sequence and also make the plots with best taxonomic classification possible. Additionally, the names of taxonomic ranks are corrected using this function.</p>
<div class="sourceCode"><pre class="sourceCode r"><code class="sourceCode r">p0.f &lt;-<span class="st"> </span><span class="kw">format_phyloseq</span>(p0)

<span class="co">#Check the taxonomy again with the formatted phyloseq object.</span>
<span class="kw">kable</span>(<span class="kw">head</span>(<span class="kw">tax_table</span>(p0.f)))</code></pre></div>
<table>
<thead>
<tr class="header">
<th></th>
<th align="left">Domain</th>
<th align="left">Phylum</th>
<th align="left">Class</th>
<th align="left">Order</th>
<th align="left">Family</th>
<th align="left">Genus</th>
<th align="left">Species</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>577110</td>
<td align="left">Bacteria</td>
<td align="left">Firmicutes</td>
<td align="left">Clostridia</td>
<td align="left">Clostridiales</td>
<td align="left">o__Clostridiales_577110</td>
<td align="left">o__Clostridiales_577110</td>
<td align="left">o__Clostridiales_577110</td>
</tr>
<tr class="even">
<td>181342</td>
<td align="left">Bacteria</td>
<td align="left">Firmicutes</td>
<td align="left">Clostridia</td>
<td align="left">Clostridiales</td>
<td align="left">Clostridiaceae</td>
<td align="left">02d06</td>
<td align="left">f__02d06_181342</td>
</tr>
<tr class="odd">
<td>581609</td>
<td align="left">Bacteria</td>
<td align="left">Firmicutes</td>
<td align="left">Clostridia</td>
<td align="left">Clostridiales</td>
<td align="left">Ruminococcaceae</td>
<td align="left">f__Ruminococcaceae_581609</td>
<td align="left">f__Ruminococcaceae_581609</td>
</tr>
<tr class="even">
<td>4341234</td>
<td align="left">Bacteria</td>
<td align="left">Firmicutes</td>
<td align="left">Clostridia</td>
<td align="left">Clostridiales</td>
<td align="left">Peptococcaceae</td>
<td align="left">Desulfotomaculum</td>
<td align="left">f__Desulfotomaculum_4341234</td>
</tr>
<tr class="odd">
<td>181348</td>
<td align="left">Bacteria</td>
<td align="left">Firmicutes</td>
<td align="left">Clostridia</td>
<td align="left">Clostridiales</td>
<td align="left">Lachnospiraceae</td>
<td align="left">Coprococcus</td>
<td align="left">f__Coprococcus_181348</td>
</tr>
<tr class="even">
<td>4467992</td>
<td align="left">Bacteria</td>
<td align="left">Firmicutes</td>
<td align="left">Bacilli</td>
<td align="left">Lactobacillales</td>
<td align="left">Streptococcaceae</td>
<td align="left">Streptococcus</td>
<td align="left">f__Streptococcus_4467992</td>
</tr>
</tbody>
</table>
<p>Check <a href="Preprocessing.md">Preprocessing</a> for additional utilitiy functions.</p>
</div>
</div>



<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

</body>
</html>