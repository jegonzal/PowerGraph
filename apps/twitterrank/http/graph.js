var metric_str = "/metrics_aggregate.json?name=Token%20Samples&rate=1&rounding=1&tlast=300"
var n = 61;
var margin = {top: 10, right: 10, bottom: 20, left: 40},
    width = 480 - margin.left - margin.right,
    height = 160 - margin.top - margin.bottom;
var x = d3.scale.linear()
  .domain([0, n - 1])
  .range([0, width]);

var y = d3.scale.linear()
  .domain([3e6, 5e6])
  .range([height, 0]);

var line = d3.svg.line()
  .x(function(d, i) { return x(i); })
  .y(function(d, i) { return y(d); });

var path = null;
var path2 = null;
var svg = null;
var matlabupdates = [];
for (var i = 0; i < n; i++)
  matlabupdates.push(2.44e6);

function get_update_rate() {
  $.getJSON(domain_str + metric_str, update_metric_graph)
      .error(function () {
        console.log("Unable to access " + domain_str + " will try again.");
      })
  .complete(function () {
    setTimeout(get_update_rate, update_interval);
  });
}

function setup_graph () {
  svg = d3.select("#update_metric").insert("svg", ":first-child")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  svg.append("defs").append("clipPath")
  .attr("id", "clip")
  .append("rect")
  .attr("width", width)
  .attr("height", height);

  svg.append("g")
  .attr("class", "x axis")
  .attr("transform", "translate(0," + height + ")")
  .call(d3.svg.axis().scale(x).orient("bottom").tickValues([]));

  svg.append("g")
  .attr("class", "y axis")
  .call(d3.svg.axis().scale(y).orient("left").tickFormat(d3.format(",s")));

  path = svg.append("g")
  .attr("clip-path", "url(#clip)")
  .append("path")

  path2 = svg.append("g")
  .attr("clip-path", "url(#clip)")
  .append("path")
}

function drawmatlab() {
  path2
  .data([matlabupdates])
  .attr("class", "line")
  .attr("d", line)
}


  function update_metric_graph(data) {
    var updates = [];
    for (var i = 0; i < data[0].record.length; i++) {
      updates.push(data[0].record[i][1]);
    }
    var ymax = Math.max(2.44e6, Math.max.apply(null, updates))  * 1.4
    var ymin = Math.min(2.44e6, Math.min.apply(null, updates))  * 0.8

    // console.log(updates);

    y = d3.scale.linear()
        .domain([ymin, ymax])
        .range([height, 0]);

    svg.select(".y.axis").call
        (d3.svg.axis().scale(y).orient("left").tickFormat(d3.format(",s")));

    drawmatlab();

      // redraw the line, and slide it to the left
      if (updates.length >= n)  {
        path.data([updates]).attr("class", "line2")
            .attr("d", line)
            .attr("transform", null)
            .transition()
            .ease("linear")
            .attr("transform", "translate(" + x(-1) + ")")
      } else {
        path.data([updates]).attr("class", "line2")
            .attr("d", line)
            .attr("transform", null)
            .transition()
            .ease("linear")
            .attr("transform", "translate(" + x(0) + ")")
      }
  }
$(function() {setup_graph();});
$(function() {get_update_rate();});
