#include <map>
#include <sstream>
#include <algorithm>
#include "paint.h"
#include "common.h"

using namespace std;


//////////////////// Paint::List ////////////////////

template<typename T>
vector<T> Paint::List<T>::to_vector() const {
	vector<T>	v;
	const List<T>	*ptr_list = this;
	while(ptr_list != NULL) {
		v.push_back(ptr_list->last);
		ptr_list = ptr_list->first.get();
	}
	
	std::reverse(v.begin(), v.end());
	return v;
}

template<typename T>
typename Paint::List<T>::ptrList Paint::List<T>::add(const ptrList ptr, T val) {
	return ptrList(new List<T>(ptr, val));
}


//////////////////// Paint::Region ////////////////////

string Paint::Region::str() const {
	stringstream	ss;
	ss << '[' << first << ", " << last << ']';
	return ss.str();
}

vector<const Paint::Region *> Paint::Region::create(const string& seq,
													const vector<double>& cMs) {
	vector<const Paint::Region *>	regions;
	size_t	last = 0U;
	Color	prev_color = -1;
	
	for(auto p = seq.begin(); p != seq.end(); ++p) {
		const Color	color = (Color)(*p) - '0';
		if(prev_color != -1 && color != prev_color) {
			const size_t	j = p - seq.begin();
			const double	cM = cMs[j-1] - cMs[last];
			regions.push_back(new Region(prev_color, last, j, cM));
			last = p - seq.begin();
		}
		prev_color = color;
	}
	const double	cM = cMs.back() - cMs[last];
	regions.push_back(new Region(prev_color, last, seq.size(), cM));
	return regions;
}


//////////////////// Paint::State ////////////////////

Paint::State::State(const std::vector<const Region *>& rs) :
						regions(rs),
						color_list(new List<Color>(regions[0]->get_color())),
						size(1U), num_continuous(1U),
						continuous_length(regions[0]->get_cM()),
						painted_len(0.0) { }

bool Paint::State::is_changable(double MIN_CROSSOVER) const {
	if(continuous_length >= MIN_CROSSOVER)
		return true;
	
	// it is painted the same color till now
	ptrColors	colors = color_list;
	while(true) {
		if(colors->back() != last_color())
			return false;
		if(colors->is_null())
			break;
		colors = colors->prev();
	}
	return true;
}

Paint::State *Paint::State::add(double cM) const {
	const bool	painted = last_color() != next_color();
	ptrColors	colors = List<Color>::add(color_list, last_color());
	const double	new_painted_length = painted_len + (painted ? cM : 0.0);
	return new State(regions, colors, size + 1, num_continuous + 1,
							continuous_length + cM, new_painted_length);
}

Paint::State *Paint::State::change(double cM) const {
	const Color	color = last_color() == 0 ? 1 : 0;
	const bool	painted = color != next_color();
	ptrColors	colors = List<Color>::add(color_list, color);
	const double	new_painted_length = painted_len + (painted ? cM : 0.0);
	return new State(regions, colors, size + 1, 1, cM, new_painted_length);
}

vector<Paint::Color> Paint::State::get_colors() const {
	vector<Color>	colors;
	vector<Color>	region_colors = color_list.get()->to_vector();
	for(auto p = region_colors.begin(); p != region_colors.end(); ++p) {
		const size_t	k = p - region_colors.begin();
		const Region	*region = *(regions.begin() + k);
		for(size_t i = 0U; i < region->get_size(); ++i)
			colors.push_back(*p);
	}
	return colors;
}

vector<const Paint::State *> Paint::State::nexts(double cM,
												double MIN_CROSSOVER) const {
	vector<const State *>	states;
	if(is_extendable())
		states.push_back(add(cM));
	if(is_changable(MIN_CROSSOVER))
		states.push_back(change(cM));
	return states;
}

Paint::ptrState Paint::State::min(const vector<ptrState>& states) {
	ptrState	min_state = states.front();
	for(auto p = states.begin() + 1; p != states.end(); ++p) {
		const State	*state = p->get();
		if(state->painted_length() < min_state.get()->painted_length())
			min_state = *p;
	}
	return min_state;
}

// minimize states
vector<Paint::ptrState> Paint::State::select(
								const vector<const State *>& states) {
	typedef pair<int,Color>	Key;
	map<Key,vector<ptrState>>	dic;
	for(auto p = states.begin(); p != states.end(); ++p) {
		Key	key((*p)->get_num_continuous(),(*p)->last_color());
		dic[key].push_back(ptrState(*p));
	}
	
	vector<ptrState>	selected_states;
	for(auto p = dic.begin(); p != dic.end(); ++p)
		selected_states.push_back(min(p->second));
	return selected_states;
}

vector<Paint::ptrState> Paint::State::next(vector<ptrState>& states,
											double cM, double MIN_CROSSOVER) {
	vector<const State *>	new_states;
	for(auto p = states.begin(); p != states.end(); ++p) {
		const auto	next_states = p->get()->nexts(cM, MIN_CROSSOVER);
		new_states.insert(new_states.end(),
							next_states.begin(), next_states.end());
	}
	return select(new_states);
}

string Paint::State::str() const {
	stringstream	ss;
	for(auto p = regions.begin(); p != regions.end(); ++p)
		ss << (*p)->str();
	ss << "\n";
	
	const vector<Color>	colors = this->get_colors();
	for(auto p = colors.begin(); p != colors.end(); ++p)
		ss << *p;
	
	return ss.str();
}


//////////////////// Paint ////////////////////

string Paint::paint(const string& seq, const vector<double>& cMs,
												double MIN_CROSSOVER) {
	const auto	regions = Region::create(seq, cMs);
	if(regions.size() <= 2U)
		return seq;
	
	auto	states = vector<ptrState>(1U, ptrState(new State(regions)));
	for(auto p = regions.begin() + 1; p != regions.end(); ++p)
		states = State::next(states, (*p)->get_cM(), MIN_CROSSOVER);
	
	ptrState	state = State::min(states);
	const auto	colors = state.get()->get_colors();
	
	stringstream	ss;
	for(auto p = colors.begin(); p != colors.end(); ++p)
		ss << *p;
	Common::delete_all(regions);
	return ss.str();
}
