/**
 * CAT12 Smart Search - AI-style FAQ search powered by Fuse.js
 * Indexes all sections of the page and provides fuzzy search with
 * relevant answer snippets. Works entirely client-side (no API needed).
 */
(function () {
  'use strict';

  // ── Configuration ──────────────────────────────────────────────
  const FUSE_OPTIONS = {
    keys: [
      { name: 'title',   weight: 0.4 },
      { name: 'content', weight: 0.6 }
    ],
    threshold: 0.35,        // 0 = exact, 1 = match anything
    ignoreLocation: true,   // search the entire string
    includeScore: true,
    includeMatches: true,
    minMatchCharLength: 3,
    findAllMatches: true
  };

  const MAX_RESULTS    = 8;
  const SNIPPET_LENGTH = 300; // characters shown per result

  // ── Build search index from page content ───────────────────────
  function buildIndex() {
    const entries = [];

    // 1) Index every <details> block (FAQ items)
    document.querySelectorAll('details').forEach(function (det) {
      var summary = det.querySelector('summary');
      var title   = summary ? summary.textContent.trim() : '';
      var clone   = det.cloneNode(true);
      var sum2    = clone.querySelector('summary');
      if (sum2) sum2.remove();
      var content = clone.textContent.trim();
      // find nearest heading
      var anchor  = findAnchor(det);
      entries.push({ title: title, content: content, anchor: anchor, type: 'faq' });
    });

    // 2) Index every major section (h2 / h4 headings with id)
    var headings = document.querySelectorAll('#content h2[id], #content h4[id]');
    headings.forEach(function (h) {
      var sectionText = collectSectionText(h);
      if (sectionText.length < 30) return; // skip very short ones
      entries.push({
        title: h.textContent.trim(),
        content: sectionText,
        anchor: h.id,
        type: 'section'
      });
    });

    return entries;
  }

  /** Collect text from a heading until the next heading of equal or higher rank */
  function collectSectionText(heading) {
    var level = parseInt(heading.tagName[1], 10);
    var parts = [];
    var el = heading.nextElementSibling;
    while (el) {
      if (/^H[1-6]$/.test(el.tagName)) {
        var elLevel = parseInt(el.tagName[1], 10);
        if (elLevel <= level) break;
      }
      parts.push(el.textContent.trim());
      el = el.nextElementSibling;
    }
    return parts.join(' ').replace(/\s+/g, ' ');
  }

  /** Walk up the DOM to find the nearest element with an id (for linking) */
  function findAnchor(el) {
    var node = el;
    while (node) {
      if (node.id) return node.id;
      // check previous siblings for headings with id
      var prev = node.previousElementSibling;
      while (prev) {
        if (prev.id) return prev.id;
        if (/^H[1-6]$/.test(prev.tagName) && prev.id) return prev.id;
        prev = prev.previousElementSibling;
      }
      node = node.parentElement;
    }
    return 'top';
  }

  /** Create a short snippet around the best‐matching region */
  function makeSnippet(text, query) {
    if (!text) return '';
    var lower  = text.toLowerCase();
    var qLower = query.toLowerCase().split(/\s+/);
    var bestPos = -1;
    // find position of first query word
    for (var i = 0; i < qLower.length; i++) {
      var idx = lower.indexOf(qLower[i]);
      if (idx !== -1) { bestPos = idx; break; }
    }
    if (bestPos === -1) bestPos = 0;
    var start = Math.max(0, bestPos - 60);
    var end   = Math.min(text.length, start + SNIPPET_LENGTH);
    var snippet = (start > 0 ? '…' : '') +
                  text.substring(start, end) +
                  (end < text.length ? '…' : '');
    return snippet;
  }

  /** Highlight matching words in text */
  function highlight(text, query) {
    if (!query) return text;
    var words = query.split(/\s+/).filter(function (w) { return w.length >= 3; });
    if (words.length === 0) return text;
    var pattern = new RegExp('(' + words.map(escapeRegex).join('|') + ')', 'gi');
    return text.replace(pattern, '<mark>$1</mark>');
  }

  function escapeRegex(s) {
    return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  }

  // ── UI ─────────────────────────────────────────────────────────
  function createSearchUI() {
    // Search box in sidebar
    var searchSection = document.createElement('section');
    searchSection.id = 'search-section';
    searchSection.innerHTML =
      '</header>' +
      '<div id="cat-search-box">' +
      '  <input type="text" id="cat-search-input" placeholder="Search for keywords in Manual" autocomplete="off" aria-label="Search the CAT manual">' +
      '  <button id="cat-search-btn" aria-label="Search"><i class="fas fa-search"></i></button>' +
      '</div>' +
      '<div id="cat-search-results" aria-live="polite"></div>';

    // Insert before the menu in sidebar
    var sidebar = document.querySelector('#sidebar > .inner');
    if (sidebar) {
      sidebar.insertBefore(searchSection, sidebar.firstChild);
    }
  }

  function renderResults(results, query) {
    var container = document.getElementById('cat-search-results');
    if (!results || results.length === 0) {
      container.innerHTML = '<p class="search-no-result">No matching answers found. Try rephrasing your question.</p>';
      return;
    }

    var html = '';
    results.slice(0, MAX_RESULTS).forEach(function (r) {
      var item    = r.item;
      var snippet = makeSnippet(item.content, query);
      var badge   = item.type === 'faq' ? '<span class="search-badge faq">FAQ</span>' : '<span class="search-badge section">Manual</span>';
      html += '<div class="search-result">' +
              '  <a href="#' + item.anchor + '" class="search-result-title">' + badge + ' ' + highlight(item.title, query) + '</a>' +
              '  <p class="search-result-snippet">' + highlight(snippet, query) + '</p>' +
              '</div>';
    });
    container.innerHTML = html;
  }

  // ── Initialisation ─────────────────────────────────────────────
  function init() {
    createSearchUI();

    var data  = buildIndex();
    var fuse  = new Fuse(data, FUSE_OPTIONS);
    var input = document.getElementById('cat-search-input');
    var btn   = document.getElementById('cat-search-btn');
    var resultsContainer = document.getElementById('cat-search-results');

    function doSearch() {
      var q = input.value.trim();
      if (q.length < 2) {
        resultsContainer.innerHTML = '';
        return;
      }
      var results = fuse.search(q);
      renderResults(results, q);
    }

    // Debounced live search
    var timer = null;
    input.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(doSearch, 300);
    });

    input.addEventListener('keydown', function (e) {
      if (e.key === 'Enter') { clearTimeout(timer); doSearch(); }
    });

    btn.addEventListener('click', function () {
      clearTimeout(timer);
      doSearch();
    });

    // Close results when clicking a result link
    resultsContainer.addEventListener('click', function (e) {
      if (e.target.closest('a')) {
        setTimeout(function () { resultsContainer.innerHTML = ''; }, 400);
      }
    });
  }

  // Wait for DOM
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
