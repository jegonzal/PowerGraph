var domain_str = "http://gldemo";
var page_str = "/wordclouds";
var update_interval = 5000;
var fill = d3.scale.category20();

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_words();
}

function get_top_words() {
    d3.json(domain_str+page_str, update_wordclouds);
    setTimeout(get_top_words, update_interval);
}

function update_wordclouds(data) {
  if (data == null) {
    console.log("Unable to access " + domain_str + " will try again.");
  }

  var topics = d3.select("#word_clouds").selectAll(".topics")
    .data(data.values)

  topics.enter()
    .append("svg")
    .attr("class", "topics")

  topics.each( function(d, i) {
    layout(d, this);
    console.log("update word clouds");
  });


  topics.exit().remove();

}

/*
 * Draw the word layout. Need this object.
 */
function draw(words) {
  var svg = d3.select(this)
              .attr("width", 300)
              .attr("height", 300);

  if (svg.select("g")[0][0] == null)
    svg.append("g");
      
  var g = svg.select("g")
    .attr("transform", "translate(150,150)")
    .selectAll("text")
    .data(words);

    g.enter().append("text")
    .style("font-size", function(d) { return d.size + "px"; })
    .style("font-family", "Impact")
    .style("fill", function(d, i) { return fill(i); })
    .attr("text-anchor", "middle")
    .attr("transform", function(d) {
      return "translate(" + [d.x, d.y] + ")rotate(" + d.rotate + ")";
    })
  .text(function(d) { return d.text; });
  g.exit().remove();
}

/* 
 * Compute word layout with data d
 * svg: the "this" object to call draw.
 * */ 
function layout(d, svg) {
    var draw_call_back = draw.bind(svg);
    d3.layout.cloud().size([300, 300])
        .words(d.map(function(d) {
          return {text: d[0], size: d[1]*0.001};
        }))
        .rotate(function() { return ~~(Math.random() * 2) * 90; })
        .font("Impact")
        .fontSize(function(d) { return d.size; })
        .on("end", draw_call_back)
        .start();
}
