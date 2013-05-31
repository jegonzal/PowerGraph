function update_user_field(data) {
  $("#userid").text("User: " + userid);
  var trainobj = data.filter(function(elem) { return (elem.training > 0) });
  var testobj = data.filter(function(elem) { return (elem.training == 0) });

  $("#mvtrain").html(
    trainobj.map(function(elem) {
    return "<div class=\"movie_query\">" +
           "<div class=\"movie_query_field\">"+elem.id+"</div>"+ 
           "<div class=\"movie_query_field\">"+elem.name + "</div>" +
           "<div class=\"movie_query_field\">"+elem.rating+"</div>" +
           "<div class=\"movie_query_field\">"+elem.pred+"</div></div>"
     }).join("\n") 
   );

  $("#mvtest").html(
    testobj.map(function(elem) {
    return "<div class=\"movie_query\">" +
           "<div class=\"movie_query_field\">"+elem.id+"</div>"+ 
           "<div class=\"movie_query_field\">"+elem.name + "</div>" +
           "<div class=\"movie_query_field\">"+elem.rating+"</div>" +
           "<div class=\"movie_query_field\">"+elem.pred+"</div></div>"
   }).join("\n"));
  // $("mvrecommend").text();
}

function get_user_info() {
    $.getJSON(domain_str + query_user_str, 
              {uid: userid},
              update_user_field)
      .error(function() { 
            console.log("Unable to access " + domain_str + " will try again.");
            jNotify("Unable to access " + domain_str + " will try again.");
        })
      .complete(function() {
            setTimeout(get_user_info, update_interval);
        });
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

get_user_info();
