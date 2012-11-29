// install slider for alpha, beta
function update_alpha(val) {
  $.get(domain_str+ldaparam_str, {"alpha": val},
      function(msg) {
        console.log("Update Alpha: " + msg);
      });
  console.log("request change alpha to :" + val);
}

function update_beta(val) {
  $.get(domain_str+ldaparam_str, {"beta": val},
        function(msg) {
          console.log("Update Beta: " + msg)
        });
  console.log("request change beta to :" + val);
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

$(function() {
  install_slider("Alpha", update_alpha);
});

$(function() {
  install_slider("Beta", update_beta);
});
