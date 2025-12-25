window.MathJax = {
  tex: {
    inlineMath: [["$","$"], ["\\(", "\\)"]],
    displayMath: [["$$","$$"], ["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    // ignoreHtmlClass: ".*|",
    ignoreHtmlClass:'none',
    // processHtmlClass: "arithmatex"
    processHtmlClass: 'none'
  }
};

document$.subscribe(() => { 
  MathJax.startup.output.clearCache()
  MathJax.typesetClear()
  MathJax.texReset()
  MathJax.typesetPromise()
})

// Inject CSS to allow long MathJax blocks to wrap when automatic linebreaking runs.
(function(){
  try{
    const css = `
.mjx-block, .MathJax_Display, .mjx-chtml, .mjx-svg {
  white-space: normal !important;
  overflow-wrap: anywhere !important;
  word-break: break-word !important;
}
.mjx-chtml { line-height: 1.2 !important; }
`;
    const style = document.createElement('style');
    style.type = 'text/css';
    style.appendChild(document.createTextNode(css));
    document.head.appendChild(style);
    console.debug('MathJax: injected wrapping CSS');
  } catch (e) { console.warn('MathJax: failed injecting CSS', e); }
})();

MathJax.startup.promise.then(() => {
  console.debug('MathJax: startup completed, automatic linebreaks:', (MathJax.chtml && MathJax.chtml.linebreaks) || (MathJax.svg && MathJax.svg.linebreaks));
}).catch(e => console.warn('MathJax: startup promise rejected', e));