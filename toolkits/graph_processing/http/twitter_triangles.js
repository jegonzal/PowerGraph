google.load("jquery", "1.5");
google.load("jqueryui", "1.7.2");
//google.load("visualization", "1");


var domain_str = "http://localhost:8090";
var domain_str = "";
var page_str =  "top_users.json";

// jsonp callback required
var twitter_addr = "http://api.twitter.com/1/users/show.json"

var current_results = [];
var users_remaining = 0;
var user_profiles = {};

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_users();
}


function refresh() {
    get_top_users();
    reloadStylesheets();
}

// Start the rendering of the UI
google.setOnLoadCallback(function() { 
    get_top_users();
});


function get_top_users() {
    jQuery.getJSON(domain_str + page_str, get_user_profiles).error(function() { 
        console.log("Unable to access " + domain_str + " will try again.");
    });
    // .complete(function() {
    //             setTimeout(get_top_users, update_interval);
    //         });
}

 


function get_user_profiles(data) {
    // save the original results
    current_results = data;
    // compute the union of all the _missing_ profiles
    jQuery.each(current_results, function(i, list) {
        console.log(list.name);
        jQuery.each(list.values, function(i, pair) {
            var id = pair[0];
            if(user_profiles[id] == undefined) { 
                users_remaining++;
                user_profiles[id] = { queried: false, is_set: false, profile: {} }; 
            }
        });
    });

    // Grab all _missing_ profiles
    jQuery.each(user_profiles, function(id, obj) {
        console.log(id);
        var query_str = twitter_addr + "?user_id=" + id;
        if(!user_profiles[id].queried) {
            console.log("Requesting: " + query_str);
            user_profiles[id].queried = true;
            jQuery.getJSON(query_str, function(data) {
                user_profiles[id].is_set = true;
                user_profiles[id].profile = data;
            }).error(function() { 
                console.log("Unable to access " + query_str + " will try again.");
            }).complete(function() { 
                users_remaining--;
                if(users_remaining == 0) { render_page(); }   
            });
        }
    });
} // end of get user profiles


function render_page() {
    var container = $("#results");
    container.empty();

    // compute the union of all the profiles
    jQuery.each(current_results, function(i, list) {
        console.log("Creating div for: " + list.name);
        var div_str = 
            "<div class=\"user_list\" id=\"" + list.name + "\">" +
            "<div class=\"title\">" + list.label + "</div>" +
            "<div class=\"contents\">"
        jQuery.each(list.values, function(i, pair) {
            var id = pair[0];
            var count = pair[1];
            if(user_profiles[id].is_set) {
                var profile = user_profiles[id].profile;
                div_str += "<div class=\"user\" id=\"" + profile.id_str + "\">" +
                    "<img class=\"user_image\" src=\"" + profile.profile_image_url + "\" / >" +
                    "<div class=\"user_info\">" +
                    "<div class=\"name\">" +
                    "<a href=\"http://twitter.com/#!/" + profile.screen_name + "\">" + 
                    profile.name + 
                    "</a>" +
                    "</div>" +
                    "<div class=\"value\">" + count + "</div>" +
                    "</div>" + "</div>";
            }
        });
        div_str += "</div></div>";
        container.append(div_str);
    });

}

function reloadStylesheets() {
    var queryString = '?reload=' + new Date().getTime();
    $('link[rel="stylesheet"]').each(function () {
        this.href = this.href.replace(/\?.*|$/, queryString);
    });
}







