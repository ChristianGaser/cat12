/**
 * CAT manual search, powered by Fuse.js (js/fuse.min.js).
 * Indexes all sections, FAQ answers and release notes of the page and shows
 * fuzzy matches with a short snippet. Works entirely client-side and offline.
 */
(function () {
  'use strict';

  // ── Configuration ──────────────────────────────────────────────
  var FUSE_OPTIONS = {
    keys: [
      { name: 'title',   weight: 0.4 },
      { name: 'content', weight: 0.6 }
    ],
    threshold: 0.35,        // 0 = exact, 1 = match anything
    ignoreLocation: true,   // search the entire string
    minMatchCharLength: 3
  };

  var MAX_RESULTS    = 10;
  var SNIPPET_LENGTH = 220; // characters shown per result

  // ── Build search index from page content ───────────────────────
  function headingText(h) {
    var clone = h.cloneNode(true);
    var anchor = clone.querySelector('.anchor');
    if (anchor) anchor.parentNode.removeChild(anchor);
    return clone.textContent.replace(/\s+/g, ' ').trim();
  }

  /** Nearest preceding heading of the given levels (e.g. 'H2') */
  function previousHeading(el, levels) {
    var node = el;
    while (node && node.id !== 'content') {
      var prev = node.previousElementSibling;
      while (prev) {
        if (levels.indexOf(prev.tagName) !== -1) return prev;
        prev = prev.previousElementSibling;
      }
      node = node.parentElement;
    }
    return null;
  }

  function buildIndex() {
    var entries = [];
    var count = {};

    // 1) Collapsible blocks: FAQ answers and release notes
    document.querySelectorAll('#content details').forEach(function (det) {
      var summary = det.querySelector('summary');
      var clone   = det.cloneNode(true);
      var sum2    = clone.querySelector('summary');
      if (sum2) sum2.parentNode.removeChild(sum2);

      var section = previousHeading(det, ['H2']);
      var sub     = previousHeading(det, ['H3']);
      var key     = section ? section.id : 'item';
      if (!det.id) {
        count[key] = (count[key] || 0) + 1;
        det.id = key + '_' + count[key];
      }
      var isFaq = section && section.id === 'faq';
      entries.push({
        title:   summary ? summary.textContent.replace(/\s+/g, ' ').trim() : '',
        content: clone.textContent.replace(/\s+/g, ' ').trim(),
        anchor:  det.id,
        type:    isFaq ? 'faq' : 'notes',
        path:    isFaq ? '' : [section, sub].filter(Boolean).map(headingText).join(' › ')
      });
    });

    // 2) Sections: the text below each heading up to the next heading
    document.querySelectorAll('#content h2[id], #content h3[id], #content h4[id], #content h5[id]').forEach(function (h) {
      var text = collectSectionText(h);
      if (text.length < 30) return; // skip very short ones
      var parents = [];
      if (h.tagName !== 'H2') {
        var h2 = previousHeading(h, ['H2']);
        if (h2) parents.push(h2);
        if (h.tagName === 'H4' || h.tagName === 'H5') {
          var h3 = previousHeading(h, ['H3']);
          if (h3) parents.push(h3);
        }
      }
      entries.push({
        title:   headingText(h),
        content: text,
        anchor:  h.id,
        type:    'section',
        path:    parents.map(headingText).join(' › ')
      });
    });

    return entries;
  }

  /** Collect text from a heading until the next heading of any level */
  function collectSectionText(heading) {
    var parts = [];
    var el = heading.nextElementSibling;
    while (el && !/^H[1-6]$/.test(el.tagName)) {
      if (el.tagName !== 'DETAILS') parts.push(el.textContent);
      el = el.nextElementSibling;
    }
    return parts.join(' ').replace(/\s+/g, ' ').trim();
  }

  // ── Fallback when Fuse.js could not be loaded ──────────────────
  function simpleSearch(data, query) {
    var words = query.toLowerCase().split(/\s+/).filter(Boolean);
    var hits = [];
    data.forEach(function (item) {
      var title = item.title.toLowerCase();
      var content = item.content.toLowerCase();
      var score = 0;
      for (var i = 0; i < words.length; i++) {
        var inTitle = title.indexOf(words[i]) !== -1;
        if (!inTitle && content.indexOf(words[i]) === -1) return;
        score += inTitle ? 3 : 1;
      }
      hits.push({ item: item, score: score });
    });
    return hits.sort(function (a, b) { return b.score - a.score; });
  }

  // ── Snippets and highlighting ──────────────────────────────────
  function escapeHtml(s) {
    return s.replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }

  function escapeRegex(s) {
    return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  }

  /** Regular expression for the query words (min. 3 characters) at word starts */
  function wordPattern(query, flags) {
    var words = query.split(/\s+/).filter(function (w) { return w.length >= 3; });
    if (words.length === 0) return null;
    return new RegExp(words.map(function (w) {
      return (/^\w/.test(w) ? '\\b' : '') + escapeRegex(w);
    }).join('|'), flags);
  }

  /** Create a short snippet around the first matching word */
  function makeSnippet(text, query) {
    if (!text) return '';
    var pattern = wordPattern(query, 'i');
    var bestPos = pattern ? text.search(pattern) : -1;
    if (bestPos === -1) bestPos = text.toLowerCase().indexOf(query.toLowerCase().split(/\s+/)[0]);
    if (bestPos === -1) bestPos = 0;
    var start = Math.max(0, bestPos - 60);
    var end   = Math.min(text.length, start + SNIPPET_LENGTH);
    return (start > 0 ? '…' : '') + text.substring(start, end) + (end < text.length ? '…' : '');
  }

  /** Escape text and wrap words of the query in <mark> */
  function highlight(text, query) {
    var pattern = wordPattern(query, 'gi');
    if (!pattern) return escapeHtml(text);
    var out = '', last = 0, m;
    while ((m = pattern.exec(text)) !== null) {
      out += escapeHtml(text.slice(last, m.index)) + '<mark>' + escapeHtml(m[0]) + '</mark>';
      last = m.index + m[0].length;
    }
    return out + escapeHtml(text.slice(last));
  }

  // ── Initialisation ─────────────────────────────────────────────
  function init() {
    var input   = document.getElementById('cat-search-input');
    var results = document.getElementById('cat-search-results');
    if (!input || !results) return;

    var data     = buildIndex();
    var fuse     = typeof Fuse === 'function' ? new Fuse(data, FUSE_OPTIONS) : null;
    var selected = -1;
    var timer    = null;

    function links() {
      return results.querySelectorAll('.search-result');
    }

    function setOpen(open) {
      results.hidden = !open;
      input.setAttribute('aria-expanded', open ? 'true' : 'false');
      if (!open) input.removeAttribute('aria-activedescendant');
    }

    function select(i) {
      var items = links();
      if (selected >= 0 && items[selected]) items[selected].classList.remove('selected');
      selected = items.length ? (i + items.length) % items.length : -1;
      if (selected >= 0) {
        var item = items[selected];
        item.classList.add('selected');
        input.setAttribute('aria-activedescendant', item.id);
        if (item.scrollIntoView) item.scrollIntoView({ block: 'nearest' });
      }
    }

    function render(hits, query) {
      selected = -1;
      if (!hits.length) {
        results.innerHTML = '<p class="search-no-result">No matching sections found. Try other or fewer words.</p>';
        setOpen(true);
        return;
      }
      results.innerHTML = hits.slice(0, MAX_RESULTS).map(function (hit, i) {
        var item  = hit.item;
        var badge = item.type === 'faq' ? '<span class="search-badge faq">FAQ</span>'
                  : item.type === 'notes' ? '<span class="search-badge">Release notes</span>' : '';
        return '<a class="search-result" role="option" id="search-result-' + i + '" href="#' + item.anchor + '">' +
               '<span class="search-result-title">' + badge + highlight(item.title, query) + '</span>' +
               (item.path ? '<span class="search-result-path">' + escapeHtml(item.path) + '</span>' : '') +
               '<span class="search-result-snippet">' + highlight(makeSnippet(item.content, query), query) + '</span>' +
               '</a>';
      }).join('');
      setOpen(true);
    }

    function doSearch() {
      var q = input.value.trim();
      if (q.length < 2) {
        results.innerHTML = '';
        setOpen(false);
        return;
      }
      render(fuse ? fuse.search(q) : simpleSearch(data, q), q);
    }

    function go(link) {
      var id = link.getAttribute('href').slice(1);
      var target = document.getElementById(id);
      setOpen(false);
      input.blur();
      if (!target) return;
      if (target.tagName === 'DETAILS') target.open = true;
      if (location.hash === '#' + id) target.scrollIntoView();
      else location.hash = id;
    }

    input.addEventListener('input', function () {
      clearTimeout(timer);
      timer = setTimeout(doSearch, 150);
    });

    input.addEventListener('focus', function () {
      if (input.value.trim().length >= 2 && results.innerHTML) setOpen(true);
    });

    input.addEventListener('keydown', function (e) {
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
        e.preventDefault();
        if (results.hidden) doSearch();
        select(selected + (e.key === 'ArrowDown' ? 1 : -1));
      } else if (e.key === 'Enter') {
        e.preventDefault();
        clearTimeout(timer);
        if (results.hidden || selected < 0) doSearch();
        var items = links();
        if (items.length) go(items[Math.max(selected, 0)]);
      } else if (e.key === 'Escape') {
        if (!results.hidden) setOpen(false);
        else { input.value = ''; input.blur(); }
      }
    });

    results.addEventListener('click', function (e) {
      var link = e.target.closest('.search-result');
      if (!link) return;
      e.preventDefault();
      go(link);
    });

    // Close the results when clicking elsewhere
    document.addEventListener('click', function (e) {
      if (!e.target.closest('.search')) setOpen(false);
    });

    // "/" or Ctrl/Cmd+K focuses the search field
    document.addEventListener('keydown', function (e) {
      var t = e.target;
      var typing = t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.isContentEditable);
      if ((e.key === '/' && !typing) || (e.key === 'k' && (e.metaKey || e.ctrlKey))) {
        e.preventDefault();
        input.focus();
        input.select();
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
