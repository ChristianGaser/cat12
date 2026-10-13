/**
 * CAT manual: colour theme toggle, navigation drawer for small screens,
 * collapsible navigation groups with highlighting of the current section,
 * section permalinks and back-to-top button. No dependencies.
 */
(function () {
  'use strict';

  var root = document.documentElement;
  var THEME_KEY = 'cat-manual-theme';  // also read by the inline script in <head>
  var darkQuery = window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)') : null;

  function forEach(list, fn) {
    Array.prototype.forEach.call(list, fn);
  }

  // ── Colour theme ───────────────────────────────────────────────
  // The page follows the system theme. A click on the toggle stores the
  // opposite theme; choosing the system theme again removes the override.
  function systemTheme() {
    return darkQuery && darkQuery.matches ? 'dark' : 'light';
  }

  function activeTheme() {
    return root.getAttribute('data-theme') || systemTheme();
  }

  function initTheme() {
    var btn = document.querySelector('.theme-toggle');
    if (!btn) return;

    function updateLabel() {
      var label = 'Switch to ' + (activeTheme() === 'dark' ? 'light' : 'dark') + ' theme';
      btn.setAttribute('aria-label', label);
      btn.setAttribute('title', label);
    }

    btn.addEventListener('click', function () {
      var next = activeTheme() === 'dark' ? 'light' : 'dark';
      try {
        if (next === systemTheme()) {
          root.removeAttribute('data-theme');
          localStorage.removeItem(THEME_KEY);
        } else {
          root.setAttribute('data-theme', next);
          localStorage.setItem(THEME_KEY, next);
        }
      } catch (e) {
        root.setAttribute('data-theme', next);  // storage blocked: apply for this visit only
      }
      updateLabel();
    });

    if (darkQuery) {
      if (darkQuery.addEventListener) darkQuery.addEventListener('change', updateLabel);
      else if (darkQuery.addListener) darkQuery.addListener(updateLabel);
    }
    updateLabel();
  }

  // ── Navigation drawer (small screens) ──────────────────────────
  function initDrawer() {
    var btn = document.querySelector('.menu-btn');
    var sidebar = document.getElementById('sidebar');
    var backdrop = document.querySelector('.backdrop');
    if (!btn || !sidebar) return;

    function setOpen(open) {
      root.classList.toggle('nav-open', open);
      btn.setAttribute('aria-expanded', open ? 'true' : 'false');
      btn.setAttribute('aria-label', open ? 'Close navigation' : 'Open navigation');
    }

    btn.addEventListener('click', function () {
      setOpen(!root.classList.contains('nav-open'));
    });
    if (backdrop) backdrop.addEventListener('click', function () { setOpen(false); });
    sidebar.addEventListener('click', function (e) {
      if (e.target.closest('a')) setOpen(false);
    });
    document.addEventListener('keydown', function (e) {
      if (e.key === 'Escape' && root.classList.contains('nav-open')) {
        setOpen(false);
        btn.focus();
      }
    });
  }

  // ── Navigation groups and current section ──────────────────────
  function initNav() {
    var sidebar = document.getElementById('sidebar');
    var nav = sidebar && sidebar.querySelector('.nav');
    if (!nav) return;

    var groups = nav.querySelectorAll('.nav-group');
    var items = [];
    var current = null;

    forEach(nav.querySelectorAll('a[href^="#"]'), function (link) {
      var target = document.getElementById(link.getAttribute('href').slice(1));
      if (target) items.push({ link: link, target: target });
    });

    function setGroup(group, open) {
      group.classList.toggle('open', open);
      var toggle = group.querySelector('.nav-toggle');
      if (toggle) toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
    }

    // Groups opened or closed by hand keep their state.
    forEach(nav.querySelectorAll('.nav-toggle'), function (toggle) {
      toggle.addEventListener('click', function () {
        var group = toggle.closest('.nav-group');
        var open = !group.classList.contains('open');
        setGroup(group, open);
        group.setAttribute('data-manual', open ? 'open' : 'closed');
      });
    });

    function keepVisible(link) {
      var top = link.offsetTop;
      var bottom = top + link.offsetHeight;
      var margin = 48;
      if (top < sidebar.scrollTop + margin) {
        sidebar.scrollTop = top - margin;
      } else if (bottom > sidebar.scrollTop + sidebar.clientHeight - margin) {
        sidebar.scrollTop = bottom - sidebar.clientHeight + margin;
      }
    }

    function activate(item) {
      if (item === current) return;
      if (current) {
        current.link.classList.remove('active');
        current.link.removeAttribute('aria-current');
      }
      current = item;
      if (!item) return;

      item.link.classList.add('active');
      item.link.setAttribute('aria-current', 'location');
      var activeGroup = item.link.closest('.nav-group');
      forEach(groups, function (group) {
        var manual = group.getAttribute('data-manual');
        group.classList.toggle('has-active', group === activeGroup);
        if (group === activeGroup) {
          if (manual !== 'closed') setGroup(group, true);
        } else if (!manual) {
          setGroup(group, false);
        }
      });
      keepVisible(item.link);
    }

    // The current section is the last one whose heading has passed the top bar.
    function update() {
      var offset = (parseFloat(getComputedStyle(root).scrollPaddingTop) || 80) + 8;
      var found = null;
      var bestTop = -Infinity;
      for (var i = 0; i < items.length; i++) {
        var top = items[i].target.getBoundingClientRect().top;
        if (top <= offset && top > bestTop) {
          bestTop = top;
          found = items[i];
        }
      }
      var atBottom = window.innerHeight + window.pageYOffset >= document.documentElement.scrollHeight - 2;
      if (atBottom && items.length) found = items[items.length - 1];
      activate(found);
    }

    var ticking = false;
    function onScroll() {
      if (ticking) return;
      ticking = true;
      window.requestAnimationFrame(function () {
        ticking = false;
        update();
      });
    }

    window.addEventListener('scroll', onScroll, { passive: true });
    window.addEventListener('resize', onScroll);
    window.addEventListener('load', update);
    update();
  }

  // ── Section permalinks ─────────────────────────────────────────
  function initPermalinks() {
    forEach(document.querySelectorAll('.doc h2[id], .doc h3[id], .doc h4[id], .doc h5[id]'), function (h) {
      var a = document.createElement('a');
      a.className = 'anchor';
      a.href = '#' + h.id;
      a.setAttribute('aria-label', 'Link to this section');
      a.textContent = '#';
      h.appendChild(a);
    });
  }

  // ── Collapsible blocks ─────────────────────────────────────────
  // Open a <details> block when it (or something inside it) is linked to,
  // and print all blocks expanded.
  function initDetails() {
    // FAQ entries get their ids from search.js after the browser has already
    // tried to scroll to the fragment, so scroll explicitly on page load.
    function openTarget(scroll) {
      var id = decodeURIComponent(location.hash.slice(1));
      var el = id && document.getElementById(id);
      var block = el && el.closest('details');
      if (!block) return;
      block.open = true;
      if (scroll === true) el.scrollIntoView();
    }

    window.addEventListener('hashchange', openTarget);
    openTarget(true);

    window.addEventListener('beforeprint', function () {
      forEach(document.querySelectorAll('details:not([open])'), function (block) {
        block.setAttribute('data-print-open', '');
        block.open = true;
      });
    });
    window.addEventListener('afterprint', function () {
      forEach(document.querySelectorAll('details[data-print-open]'), function (block) {
        block.removeAttribute('data-print-open');
        block.open = false;
      });
    });
  }

  // ── MATLAB links and back-to-top button ────────────────────────
  function initMisc() {
    forEach(document.querySelectorAll('a.function'), function (a) {
      if (!a.title) a.title = 'Runs in MATLAB';
    });

    var btn = document.querySelector('.back-to-top');
    if (!btn) return;
    function update() {
      btn.classList.toggle('visible', window.pageYOffset > 900);
    }
    window.addEventListener('scroll', update, { passive: true });
    update();
  }

  function init() {
    initTheme();
    initDrawer();
    initPermalinks();
    initNav();
    initDetails();
    initMisc();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
