var domain_str = "http://bros.ml.cmu.edu:8090";
var wordcloudpage_str = "/wordclouds";
var clickpage_str= "/click";
var addtopic_str= "/addtopic";
var ldaparam_str= "/ldaparam";
var lockword_str = "/lockword";
var update_interval = 5000;
var fill = d3.scale.category20();
var basefont = 20;
var lockedwords = {};
var maxfont = 40;
var reorder = false;
var numtopics = 0;
var alpha = 0.15;
var beta = 0.15;

function lockword(topic, word) {
  $.get(domain_str+lockword_str, {"word": word, "topic": topic},
        function(msg) {
          console.log("Lock word: " + msg);
          // jSuccess("Lock \"" + word + "\"" + " for topic " + topic);
          if (msg == "ok") {
            if (lockedwords[topic] == undefined) {
              lockedwords[topic]={}
            }
            lockedwords[topic][word] = true;
            jSuccess("Locked \""+word+"\" for topic " + topic);
          } else if (msg == "Unable to find word") {
            jNotify(msg + ": " + word);
          } else if (msg == "Invalid Topic number") {
            jNotify(msg + ": " + topic);
          } else {
            jNotify(msg);
          }
        });
  console.log("request lock word:" + word + " for topic " + topic);
}

function on_topic_click(i) {
  // display input box
  console.log("click topic " + i);
  var textbox = $("#topic"+i+" input");
  textbox.show().focus();
  textbox.keyup(function(e) {
     // the enter key
    if (e.keyCode == 13) {
      var word = $(this).val()
      on_topic_submit(i,word);
      $(this).blur();
    }
    // esc
    if (e.keyCode == 27)  {
      $(this).blur();
    }
  });
}

function on_textbox_exit(e) {
  $(this).val('').hide().unbind("keyup");
}

function on_topic_submit(topicid, word) {
  lockword(topicid, word);
}

function on_add_topic () {
  var seeds = $("#seedinput").val()
  add_topics(seeds)
}

function add_topics(words) {
  $.get(domain_str+addtopic_str, {"seed": words},
        function(msg) {
          console.log("Add topic: " + msg);
          if (msg == "ok") {
            lockedwords[numtopics] = {};
            words = words.split(",");
            for (var i = 0; i < words.length; i++)
              lockedwords[numtopics][words[i]] = true;
            jSuccess("Added new topic " + (numtopics));
          } else {
            jNotify(msg);
          }
        });
  console.log("request add topic with seeds: " + words);
}


//================ LDA Visualizer ===========================
function get_top_words() {
    $.getJSON(domain_str + wordcloudpage_str, update_wordclouds)
      .error(function() { 
            console.log("Unable to access " + domain_str + " will try again.");
            jNotify("Unable to access " + domain_str + " will try again.");
        })
      .complete(function() {
            setTimeout(get_top_words, update_interval);
        });
}

function update_wordclouds(data) {
  // update other stuff
  $("#ndocs").text(data.ndocs);
  $("#ntokens").text(data.ntokens);
  $("#nwords").text(data.nwords);
  $("#ntopics").text(data.ntopics);
  if (data.alpha != alpha  || data.beta != beta) {
    $("#Alpha .slider").slider("option", "value", data.alpha * 100);
    $("#Beta .slider").slider("option", "value", data.beta * 100);
    alpha = data.alpha;
    beta = data.beta;
  }

  var topics = data.values;

  // update the word cloud for each topic
  for (var i = 0; i < topics.length; i++)
    update_topic(i, topics[i]);


  if (reorder) {
    // rearrange the topics by its popularity 
    var weights = [];
    for (var i = 0; i< topics.length; i++) {
      sum = 0;
      for (var j = 0; j < topics[i].length; j++) {
        sum += topics[i][j][1];
      }
      weights.push(sum);
    }
    d3.selectAll("#word_clouds .topic").data(weights).sort(compare);
    reorder = false;
  }
}

function update_topic(i, data) {
  var topic = d3.select("#word_clouds").select("#topic"+i);
    // .select("g");
  var minval = Math.max(mincount(data), 1);

  // Initialize wordcloud box for each topic
  if (topic.empty())  {
    topic = d3.select("#word_clouds")
              .insert("div", ":first-child")
              .attr("id", "topic"+i)
              .attr("class", "topic")
    topic.on("click", function(d, j) { 
      on_topic_click(i)
    });
    // install the hidden input box
    topic.append("input")
      .attr("type", "textbox")
      .style("size", 10);
    $("#topic"+i+" input").hide();
    $("#topic"+i+" input").focusout(on_textbox_exit);
    numtopics+=1;
    console.log("Add topic: " + (numtopics-1));
  }

  var wc = topic.selectAll("div").data(data, function(d) {return d[0];});

  // update words that already exist
  wc.attr("class", 
          function(d) {
            if (lockedwords[i] && lockedwords[i][d[0]]) {
              return "word locked";
            } else {
              return "word update";
            }
          })
    .transition()
    .duration(1000)
    .delay(200)
    .style("font-size", function(d) {return Math.min(maxfont, basefont + d[1]/minval)+"px";})

  // insert new words in the topic
    var enters = wc.enter().insert("div", ":nth-child(2)")
                   .attr("class", function(d){
                    if (lockedwords[i] && lockedwords[i][d[0]]) {
                      return "word locked";
                    } else {
                      return "word enter";
                    }
                   })
                    .text(function(d) {return d[0];})
                    .transition()
                    .duration(1000)
                    .delay(100)
                    .style("font-size", function(d) {return (Math.min(maxfont, basefont + d[1]/minval)+"px");})

  // remove words
  wc.exit().attr("class", "exit word")
    .transition()
    .duration(1000)
    .style("font-size", "1px")
    .remove();
}


// helper function to get the smallest weight 
function mincount(data) {
  var min = 1e20;
  for (var i = 0; i < data.length; i++) {
    if (min > data[i][1]) {
      min = data[i][1];
    }
  }
  return min;
}


// helper function for generic comparison
function compare(a,b) {
  if (a < b) return 1;
  if (a == b) return 0;
  if (a > b) return -1;
}

get_top_words();
