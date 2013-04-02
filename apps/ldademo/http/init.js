// install slider for alpha, beta
function update_alpha(val) {
  $.get(domain_str+ldaparam_str, {"alpha": val},
      function(msg) {
        console.log(msg);
        jSuccess(msg);
      });
  console.log("request change alpha to :" + val);
}

function update_beta(val) {
  $.get(domain_str+ldaparam_str, {"beta": val},
        function(msg) {
          console.log(msg)
          jSuccess(msg);
        });
  console.log("request change beta to :" + val);
}

function resetword() {
  $.get(domain_str+lockword_str, {reset: 0},
        function(msg) {
          console.log("Reset words: " + msg);
          if (msg == "reset") {
            lockedwords = {};
            jSuccess("Reset locked words");
          } else {
            jNotify (msg);
          }
        });
  console.log("request reset words");
}

function reorder_topics() {
  reorder = true;
  jNotify("Topics will be reodered in the next update...");
} 

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_words();
}

function install_slider(name, handler) {
  $("#"+name+" .slider").slider({
    min: 1,
    max: 100,
    value: 15,
  create: function() {
    var value = $("#"+name+" .slider").slider("option","value");
    $("#"+name).find(".label").text(name+": " + value/100);
  },
  change: function() {
    var value = $("#"+name+" .slider").slider("option","value");
    $("#"+name).find(".label").text(name+": " +value/100);
    handler(value/100);
  },
  slide: function() {
    var value = $("#"+name+" .slider").slider("option","value");
    value /= 100;
    $("#"+name).find(".label").text(name+": "+value);
  }
  });
}


function install_button(name, handler) {
  $("#"+name).button().click(handler)
}

$(function() {
  install_slider("Alpha", update_alpha);
});

$(function() {
  install_slider("Beta", update_beta);
});

$(function() {
  install_button("resetwords", resetword);
});

$(function() {
  install_button("reordertopics", reorder_topics);
});

$(function() {
  install_button("addtopic", on_add_topic); 
});
