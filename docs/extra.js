// weathergenr pkgdown overrides.
//
// pkgdown picks this file up automatically -- it needs no entry in
// _pkgdown.yml -- and copies it to docs/extra.js.
//
// Sole job: put the point after section numbers in the sidebar table of
// contents, so it matches the headings.
//
// pkgdown.js builds that TOC at run time with bootstrap-toc, which labels each
// entry using jQuery .text(). CSS generated content is not part of .text(), so
// the `.header-section-number::after` rule in extra.css reaches the headings
// but not the TOC. Inserting the point as real text here is what keeps the two
// in step.
//
// Runs on `load` rather than DOMContentLoaded so pkgdown.js has already built
// the TOC, and is written to be idempotent: the pattern requires whitespace
// straight after the number, so an entry that already carries a point is left
// alone however many times this runs.
window.addEventListener('load', function () {
  document.querySelectorAll('#toc .nav-link').forEach(function (link) {
    link.textContent = link.textContent.replace(/^(\d+(?:\.\d+)*)\s/, '$1. ');
  });
});
