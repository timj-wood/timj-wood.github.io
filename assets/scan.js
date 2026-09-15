// Draws a stylised genome scan (Manhattan-plot style) into any <svg class="scan">.
// Deterministic, so it looks the same on every load.
(function () {
  var svg = document.querySelector('svg.scan');
  if (!svg) return;
  var W = 1200, H = 150;
  svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);

  var seed = 7;
  function rand() {
    seed = (seed * 16807) % 2147483647;
    return (seed - 1) / 2147483646;
  }

  var chromCount = 18, lengths = [], total = 0;
  for (var c = 0; c < chromCount; c++) {
    var l = 0.5 + rand() * 1.2;
    lengths.push(l);
    total += l;
  }

  var peaks = [
    { x: 0.17, w: 0.02, h: 0.9 },
    { x: 0.46, w: 0.015, h: 0.7 },
    { x: 0.78, w: 0.025, h: 1.0 }
  ];

  var out = '', x0 = 0;
  for (var i = 0; i < chromCount; i++) {
    var w = (lengths[i] / total) * W;
    var fill = i % 2 === 0 ? '#6b4e8c' : '#4f6b4a';
    var n = Math.round(w * 0.9);
    for (var k = 0; k < n; k++) {
      var x = x0 + rand() * w;
      var t = x / W;
      var y = rand() * rand() * 0.28;
      for (var p = 0; p < peaks.length; p++) {
        var d = (t - peaks[p].x) / peaks[p].w;
        y += peaks[p].h * Math.exp(-d * d) * rand();
      }
      if (y > 1) y = 1;
      var cy = H - 14 - y * (H - 30);
      var r = 1.2 + y * 1.3;
      var op = 0.35 + y * 0.6;
      out += '<circle cx="' + x.toFixed(1) + '" cy="' + cy.toFixed(1) +
             '" r="' + r.toFixed(1) + '" fill="' + fill +
             '" fill-opacity="' + op.toFixed(2) + '"/>';
    }
    x0 += w;
  }
  out += '<line x1="0" y1="' + (H - 10) + '" x2="' + W + '" y2="' + (H - 10) +
         '" stroke="#cfd6cf" stroke-width="1"/>';
  svg.innerHTML = out;
})();
