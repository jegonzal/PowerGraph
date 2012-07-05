google.load("jquery", "1.5");



var domain_str = "http://localhost:8090";
var domain_str = "";
var page_str =  "top_users.json";

// jsonp callback required
var twitter_addr = "http://api.twitter.com/1/users/lookup.json"

var current_results = [];
var user_profiles = {};

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_users();
}


function refresh() {
    get_top_users();
    reloadStylesheets();
    render_page();
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
                user_profiles[id] = { queried: false, is_set: false, profile: {} }; 
            }
        });
    });

    var id_list = "";
    var id_list_len = 0;
    // Grab all _missing_ profiles
    jQuery.each(user_profiles, function(id, obj) {
        console.log(id);
        if(!user_profiles[id].queried) {
            console.log("Requesting: " + id);
            user_profiles[id].queried = true;
            id_list += id;
            id_list_len++;
            if(id_list_len >= 99) {
                jQuery.getJSON(twitter_addr + "?callback=?", {user_id: id_list}, process_ids);
                id = "";
                id_list_len = 0;
            } else { id_list += ","; }
        } 
    });
    if(id_list_len > 0) {
        jQuery.getJSON(twitter_addr + "?callback=?", {user_id: id_list}, process_ids);
        id = "";
        id_list_len = 0;
    }
} // end of get user profiles


function process_ids(data) {
    jQuery.each(data, function(i, profile) {
        var id = profile.id;
        user_profiles[id].is_set = true;
        user_profiles[id].profile = profile;
    });;
    render_page();
}

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







