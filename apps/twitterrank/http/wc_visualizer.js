var domain_str = "http://bros.ml.cmu.edu:8090";
var wordcloudpage_str = "/wordclouds";
var clickpage_str= "/click";
var ldaparam_str= "/ldaparam";
var lockword_str = "/lockword";
var update_interval = 5000;
var fill = d3.scale.category20();
var basefont = 20;

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_words();
}

function lockword(topic, word) {
  $.get(domain_str+lockword_str, {"word": word, "topic": topic},
        function(msg) {
          console.log("Lock word: " + msg);
          // jSuccess("Lock \"" + word + "\"" + " for topic " + topic);
          if (msg == "ok") {
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
  // jNotify("Request locking \"" + word + "\"" + " for topic " + topic);
}

function resetword() {
  $.get(domain_str+lockword_str, {reset: 0},
        function(msg) {
          console.log("Reset words: " + msg);
          if (msg == "reset") {
            jSuccess("Reset word");
          } else {
            jNotify (msg);
          }
        });
  console.log("request reset words");
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
  // clean and hide input box
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

  var topics = data.values;
  for (var i = 0; i < topics.length; i++)
    update_topic(i, topics[i]);
}



function update_topic(i, data) {
  var topic = d3.select("#word_clouds").select("#topic"+i);
    // .select("g");
  var minval = Math.max(mincount(data), 1);

  // Initialize wordcloud box for each topic
  if (topic.empty())  {
    topic = d3.select("#word_clouds")
              .append("div")
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
  }

  var wc = topic.selectAll("div").data(data, function(d) {return d[0];});

  // update words that already exist
  wc.attr("class", "word update")
    .transition()
    .duration(1000)
    .delay(200)
    .style("font-size", function(d) {return (basefont + d[1]/minval)+"px";})

  // insert new words in the topic
  var enters = wc.enter().insert("div", ":nth-child(2)").attr("class", "word enter")
  .text(function(d) {return d[0];})
  .transition()
  .duration(1000)
  .delay(100)
  .style("font-size", function(d) {return (basefont + d[1]/minval)+"px";})


  // Sort words by size
  // wc.sort(function(a, b) {
  //   if (a[1] == b[1])
  //     return a[0] < b[0];
  //   return a[1] < b[1];
  // });

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

get_top_words();
