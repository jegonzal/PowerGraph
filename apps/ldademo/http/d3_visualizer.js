var domain_str = "http://bros.ml.cmu.edu:8090";
var wordcloudpage_str = "/wordclouds";
var pagerankpage_str= "/pagerank";
var clickpage_str= "/click";
var update_interval = 5000;
var fill = d3.scale.category20();
var weights = [];

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_words();
    get_top_pages();
}

//================ LDA Visualizer ===========================
function get_top_words() {
    d3.json(domain_str+wordcloudpage_str, update_wordclouds);
    setTimeout(get_top_words, update_interval);
}

function sort_topics(data) {
  // augement the data 
  var data2 = [];
  for (var i = 0; i < data.length; i++) {
    // id, weight, topics
    data2[i] = {"id":i, "weight":weights[i], "words":data[i]};
  }

  data2.sort( function (a, b) { 
    if(a.weight == b.weight) return 0;
    return a.weight > b.weight; });

  console.log("topic ordering: " + data2.map(function(x) {return x["id"];}));
  return data2;
}

function update_wordclouds(data) {
  if (data == null) {
    console.log("Unable to access " + domain_str + " will try again.");
    return;
  }

  while (weights.length < data.values.length) weights.push(0);

  var topics = d3.select("#word_clouds").selectAll(".topics")
    .data(sort_topics(data.values))

  topics.enter()
    .append("svg")
    .attr("class", "topics")
    .on("click", function(d, i) { 
      $.get(domain_str+clickpage_str, {topic: d.id});
      alert ("Boost preference to topic" + d.id);
    })

  topics.each( function(d, i) {
    layout(d["words"], this);
    console.log("update word clouds");
  });

  topics.exit().remove();
}

/* 
 * Compute word layout with data d = ['river':1, 'mountain':2, ...]
 * svg: the "this" object to call draw.
 * */ 
function layout(d, svg) {
    var draw_call_back = draw.bind(svg);
    d3.layout.cloud().size([300, 300])
        .words(d.map(function(d) {
          return {text: d[0], size: 10+Math.random()*90};
        }))
        .rotate(function() { return ~~(Math.random() * 2) * 90; })
        .font("Impact")
        .fontSize(function(d) { return d.size; })
        .on("end", draw_call_back)
        .start();
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



//================ Pagerank Visualizer ===========================
function get_top_pages() {
    d3.json(domain_str+pagerankpage_str, update_top_pages);
    setTimeout(get_top_pages, update_interval);
}

function update_top_pages(data) {
  if (data == null) {
    console.log("Unable to access " + domain_str + " will try again.");
    return;
  }
  console.log("update pagerank vector");

  d = [Math.random(), Math.random(), Math.random()];
  var toppages= d3.select("#pagerank").selectAll("p")
    .data(data.values, function(d) { return d;});

  toppages.enter().append("p")
    .text(function(d) { 
      // vid: rank [numdocs: topic1, topic2, ... ]
      return "id:" + d[0] + " rank:" + d[1] + " topics(" + d[3] +")";
    });
  toppages.exit().remove();

  // update the global weights for each topic
  for (var j = 0; j < weights.length; j++) {
    for (var i = 0; i < data.values.length; i++) {
      weights[j] += data.values[i][3][j];
    }
  }
  console.log("topic weights: " + weights);
}
